/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *   IP, Inti Pelupessy, Netherlands eScience Center, i.pelupessy@esciencecenter.nl
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include "RooFit/MultiProcess/worker.h"

#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/types.h"
#include "RooFit/MultiProcess/Messenger.h"
#include "RooFit/MultiProcess/Job.h"
#include "RooFit/MultiProcess/util.h"
#include "RooFit/MultiProcess/ProcessTimer.h"
#include "RooFit/MultiProcess/Config.h"

#include <string>
#include <unistd.h> // getpid, pid_t
#include <cerrno>   // EINTR
#include <csignal>  // sigprocmask etc

namespace RooFit {
namespace MultiProcess {

static bool worker_loop_running = false;

bool is_worker_loop_running()
{
   return worker_loop_running;
}

/// \brief The worker processes' event loop
///
/// Asks the queue process for tasks, polls for incoming messages from other
/// processes and handles them.
void worker_loop()
{
   assert(JobManager::instance()->process_manager().is_worker());
   worker_loop_running = true;
   Q2W message_q2w;

   // use a flag to not ask twice
   bool dequeue_acknowledged = true;

   ZeroMQPoller poller;
   std::size_t mw_sub_index;

   std::tie(poller, mw_sub_index) = JobManager::instance()->messenger().create_worker_poller();

   // Before blocking SIGTERM, set the signal handler, so we can also check after blocking whether a signal occurred
   // In our case, we already set it in the ProcessManager after forking to the queue and worker processes.

   sigset_t sigmask;
   sigemptyset(&sigmask);
   sigaddset(&sigmask, SIGTERM);
   sigprocmask(SIG_BLOCK, &sigmask, &JobManager::instance()->messenger().ppoll_sigmask);

   // Before doing anything, check whether we have received a terminate signal while blocking signals!
   // In this case, we also do that in the while condition.
   while (!ProcessManager::sigterm_received()) {
      try { // watch for error from ppoll (which is called inside receive functions) caused by SIGTERM from master

         // try to dequeue a task
         if (dequeue_acknowledged) { // don't ask twice
            JobManager::instance()->messenger().send_from_worker_to_queue(W2Q::dequeue);
            dequeue_acknowledged = false;
         }

         std::chrono::time_point<std::chrono::steady_clock> polling_begin = std::chrono::steady_clock::now();
         // wait for handshake from queue or update from SUB socket
         if (RooFit::MultiProcess::Config::getTimingAnalysis()) ProcessTimer::start_timer("worker:polling");
         auto poll_result = poller.ppoll(-1, &JobManager::instance()->messenger().ppoll_sigmask);
         if (RooFit::MultiProcess::Config::getTimingAnalysis()) ProcessTimer::end_timer("worker:polling");

         //if (MultiProcess::Config::isInLinesearch_) std::cout << getpid() << " polling wait: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - polling_begin).count();
         //if (MultiProcess::Config::isInLinesearch_) std::cout << " | current time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - MultiProcess::Config::begin).count() << std::endl;	
         // because the poller may now have a waiting update from master over the SUB socket,
         // but the queue socket could be first in the poll_result vector, and during handling
         // of a new task it is possible we need to already receive the updated state over SUB,
         // we have to then flip this boolean so that in the for loop when we reach the SUB
         // socket's result, we can skip it (otherwise we will hang there, because no more
         // updated state will be coming):
         bool skip_sub = false;
         // then process incoming messages from sockets
         for (auto readable_socket : poll_result) {
            // message comes from the master-worker SUB socket (first element):
            if (readable_socket.first == mw_sub_index) {
               if (!skip_sub) {
                  std::chrono::time_point<std::chrono::steady_clock> worker_update_state_begin = std::chrono::steady_clock::now();
                  auto job_id = JobManager::instance()->messenger().receive_from_master_on_worker<std::size_t>();
                  JobManager::get_job_object(job_id)->update_state();
                  if (MultiProcess::Config::isInLinesearch_) std::cout << getpid() << " WORKER update state: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_begin).count();
                  if (MultiProcess::Config::isInLinesearch_) std::cout << " | current time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - MultiProcess::Config::begin).count() << std::endl;	
               }
            } else { // from queue socket
               message_q2w = JobManager::instance()->messenger().receive_from_queue_on_worker<Q2W>();
               switch (message_q2w) {
               case Q2W::dequeue_rejected: {
                  dequeue_acknowledged = true;
                  break;
               }
               case Q2W::dequeue_accepted: {
                  dequeue_acknowledged = true;
                  std::chrono::time_point<std::chrono::steady_clock> worker_receiving_begin = std::chrono::steady_clock::now();
                  auto job_id = JobManager::instance()->messenger().receive_from_queue_on_worker<std::size_t>();
                  auto state_id = JobManager::instance()->messenger().receive_from_queue_on_worker<State>();
                  auto task_id = JobManager::instance()->messenger().receive_from_queue_on_worker<Task>();
                  if (MultiProcess::Config::isInLinesearch_) std::cout << getpid() << " WORKER receiving: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_receiving_begin).count();
                  if (MultiProcess::Config::isInLinesearch_) std::cout << " | current time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - MultiProcess::Config::begin).count() << std::endl;	


                  // while loop, because multiple jobs may have updated state coming
                  while (state_id != JobManager::get_job_object(job_id)->get_state_id()) {
                     skip_sub = true;
                     auto job_id_for_state =
                        JobManager::instance()->messenger().receive_from_master_on_worker<std::size_t>();
                     JobManager::get_job_object(job_id_for_state)->update_state();
                  }

                  std::chrono::time_point<std::chrono::steady_clock> eval_task_begin = std::chrono::steady_clock::now();		

                  if (RooFit::MultiProcess::Config::getTimingAnalysis()) ProcessTimer::start_timer("worker:eval_task:" + std::to_string(task_id));
                  JobManager::get_job_object(job_id)->evaluate_task(task_id);
                  if (RooFit::MultiProcess::Config::getTimingAnalysis()) ProcessTimer::end_timer("worker:eval_task:" + std::to_string(task_id));
                  JobManager::get_job_object(job_id)->send_back_task_result_from_worker(task_id);

                  if (MultiProcess::Config::isInLinesearch_) std::cout << getpid() << " WORKER evaluate task ";		
                  if (MultiProcess::Config::isInLinesearch_) std::cout << "(wtime): " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - eval_task_begin).count();
                  if (MultiProcess::Config::isInLinesearch_) std::cout << " | current time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - MultiProcess::Config::begin).count() << std::endl;

                  break;
               }
               }
            }
         }

      } catch (ZMQ::ppoll_error_t &e) {
         zmq_ppoll_error_response response;
         try {
            response = handle_zmq_ppoll_error(e);
         } catch (std::logic_error &) {
            printf("worker loop at PID %d got unhandleable ZMQ::ppoll_error_t\n", getpid());
            throw;
         }
         if (response == zmq_ppoll_error_response::abort) {
            break;
         } else if (response == zmq_ppoll_error_response::unknown_eintr) {
            printf("EINTR in worker loop at PID %d but no SIGTERM received, continuing\n", getpid());
            continue;
         } else if (response == zmq_ppoll_error_response::retry) {
            printf("EAGAIN from ppoll in worker loop at PID %d, continuing\n", getpid());
            continue;
         }
      } catch (zmq::error_t &e) {
         printf("unhandled zmq::error_t (not a ppoll_error_t) in worker loop at PID %d with errno %d: %s\n", getpid(),
                e.num(), e.what());
         throw;
      }
   }

   if (RooFit::MultiProcess::Config::getTimingAnalysis()) ProcessTimer::write_file();
   
   // clean up signal management modifications
   sigprocmask(SIG_SETMASK, &JobManager::instance()->messenger().ppoll_sigmask, nullptr);

   worker_loop_running = false;
}

} // namespace MultiProcess
} // namespace RooFit
