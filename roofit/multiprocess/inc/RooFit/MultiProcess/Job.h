/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *   IP, Inti Pelupessy, Netherlands eScience Center, i.pelupessy@esciencecenter.nl
 *
 * Copyright (c) 2016-2019, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_MultiProcess_Job_decl
#define ROOT_ROOFIT_MultiProcess_Job_decl

#include <string>
#include <zmq.hpp>

namespace RooFit {
namespace MultiProcess {

// forward declaration
class JobManager;

/*
 * @brief interface class for defining the actual work that must be done
 *
 * Think of "job" as in "employment", e.g. the job of a baker, which
 * involves *tasks* like baking and selling bread. The Job must define the
 * tasks through its execution (evaluate_task), based on a task index argument.
 *
 * Classes inheriting from Job must implement the pure virtual methods:
 * - void evaluate_task(std::size_t task)
 * - void send_back_task_result_from_worker(std::size_t task)
 * - void receive_task_result_on_master(const zmq::message_t & message)
 *
 * An example/reference implementation can be found in test_Job.cxx.
 *
 * Most Jobs will also want to override the virtual update_state() function.
 * This function can be used to send and receive state from master to worker.
 * In the worker loop, when something is received over the ZeroMQ "SUB" socket,
 * update_state() is called to put the received data into the right places,
 * thus updating for instance parameter values on the worker that were updated
 * since the last call on the master side.
 *
 * ## Message protocol
 *
 * One simple rule must be upheld for the messages that the implementer will
 * send with 'send_back_task_result_from_worker' and 'update_state': the first
 * part of the message must always be the 'Job''s ID, stored in 'Job::id'.
 * The rest of the message, i.e. the actual data to be sent, is completely up
 * to the implementation. Note that on the receiving end, i.e. in the
 * implementation of 'receive_task_result_on_master', one will get the whole
 * message, but the 'Job' ID part will already have been identified in the
 * 'JobManager', so one needn't worry about it further inside
 * 'Job::receive_task_result_on_master' (it is already routed to the correct
 * 'Job'). The same goes for the receiving end of 'update_state', except that
 * update_state is routed from the 'worker_loop', not the 'JobManager'.
 *
 * ## Implementers notes
 *
 * The type of result from each task is strongly dependent on the Job at hand
 * and so Job does not provide a default results member. It is up to the
 * inheriting class to implement this in the above functions. We would have
 * liked a template parameter task_result_t, so that we could also provide a
 * default "boilerplate" calculate function to show a typical Job use-case of
 * all the above infrastructure. This is not trivial, because the JobManager
 * has to keep a list of Job pointers, so if there would be different template
 * instantiations of Jobs, this would complicate this list.
 *
 * Child classes should refrain from direct access to the JobManager instance
 * (through JobManager::instance), but rather use the here provided
 * Job::get_manager(). This function starts the worker_loop on the worker when
 * first called, meaning that the workers will not be running before they
 * are needed.
 */
class Job {
public:
   explicit Job();
   Job(const Job &other);

   ~Job();

   virtual void evaluate_task(std::size_t task) = 0;
   virtual void update_state();

   virtual void send_back_task_result_from_worker(std::size_t task) = 0;
   virtual bool receive_task_result_on_master(const zmq::message_t & message) = 0;

   void gather_worker_results();

protected:
   JobManager *get_manager();

   std::size_t id;

private:
   // do not use _manager directly, it must first be initialized! use get_manager()
   JobManager *_manager = nullptr;
};

} // namespace MultiProcess
} // namespace RooFit

#endif // ROOT_ROOFIT_MultiProcess_Job_decl
