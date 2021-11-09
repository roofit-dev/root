/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *   IP, Inti Pelupessy, Netherlands eScience Center, i.pelupessy@esciencecenter.nl
 *
 * Copyright (c) 2016-2021, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */

#include "RooFit/MultiProcess/Messenger.h"
#include "RooFit/MultiProcess/util.h"

#include <csignal> // sigprocmask etc
#include <cstdio>  // sprintf

namespace RooFit {
namespace MultiProcess {

void set_socket_immediate(ZmqLingeringSocketPtr<> &socket)
{
   int optval = 1;
   socket->set(zmq::sockopt::immediate, optval);
}

/** \class Messenger
 *
 * \brief Manages ZeroMQ sockets and wraps send and receive calls
 *
 * This class is used for all interprocess communication between the master,
 * queue and worker processes. It sets up ZeroMQ sockets between all processes
 * over IPC socket files stored in /tmp on the filesystem.
 *
 * Several sockets are used for communication between different places for
 * different purposes:
 * - Master and queue processes each have a PUSH-PULL socket pair to directly
 *   send/receive data between only the master and queue processes. This is
 *   currently used mainly for sending tasks to the queue from master. The
 *   socket from queue back to master is used only to test connections and may
 *   be removed in the future.
 * - The queue process also has a PUSH-PULL socket pair with each worker
 *   process. These are used by the workers to obtain tasks from the queue.
 * - The master has a PUB socket that the workers subscribe to with SUB
 *   sockets. These are used to update state. Note that to ensure robust
 *   reception of all messages on the SUB socket, it's important to send over
 *   state in as little messages as possible. For instance, it's best to send
 *   arrays over in a single big message instead of sending over each element
 *   separately. This also improves performance, since each message has some
 *   fixed overhead.
 * - Each worker has a PUSH socket connected to a PULL socket on master that
 *   is used to send back task results from workers to master in
 *   'JobManager::retrieve()'.
 *
 * @param process_manager ProcessManager instance which manages the master,
 *                        queue and worker processes that we want to set up
 *                        communication for in this Messenger.
 */

Messenger::Messenger(const ProcessManager &process_manager)
{
   sigemptyset(&ppoll_sigmask);

   // high water mark for master-queue sending, which can be quite a busy channel, especially at the start of a run
   int hwm = 0;
   // create zmq connections (zmq context is automatically created in the ZeroMQSvc class and maintained as singleton)
   // and pollers where necessary
   try {
      if (process_manager.is_master()) {
         mq_push_.reset(zmqSvc().socket_ptr(zmq::socket_type::push));
         auto rc = zmq_setsockopt(*mq_push_, ZMQ_SNDHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         mq_push_->bind("ipc:///tmp/roofitMP_from_master_to_queue");

         mq_push_poller_.register_socket(*mq_push_, zmq::event_flags::pollout);

         mq_pull_.reset(zmqSvc().socket_ptr(zmq::socket_type::pull));
         rc = zmq_setsockopt(*mq_pull_, ZMQ_RCVHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         mq_pull_->bind("ipc:///tmp/roofitMP_from_queue_to_master");

         mq_pull_poller_.register_socket(*mq_pull_, zmq::event_flags::pollin);

         mw_pub_.reset(zmqSvc().socket_ptr(zmq::socket_type::pub));
         rc = zmq_setsockopt(*mw_pub_, ZMQ_SNDHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         mw_pub_->bind("ipc:///tmp/roofitMP_from_master_to_workers");

         wm_pull_.reset(zmqSvc().socket_ptr(zmq::socket_type::pull));
         rc = zmq_setsockopt(*wm_pull_, ZMQ_RCVHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         wm_pull_->bind("ipc:///tmp/roofitMP_from_workers_to_master");
         wm_pull_poller_.register_socket(*wm_pull_, zmq::event_flags::pollin);

         close_MQ_on_destruct_ = true;

         // make sure all subscribers are connected
         ZmqLingeringSocketPtr<> subscriber_ping_socket{zmqSvc().socket_ptr(zmq::socket_type::rep)};
         subscriber_ping_socket->bind("ipc:///tmp/roofitMP_subscriber_ping_socket");
         ZeroMQPoller subscriber_ping_poller;
         subscriber_ping_poller.register_socket(*subscriber_ping_socket, zmq::event_flags::pollin);
         std::size_t N_subscribers_confirmed = 0;
         while (N_subscribers_confirmed < process_manager.N_workers()) {
            zmqSvc().send(*mw_pub_, false);
            auto poll_results = subscriber_ping_poller.poll(0);
            for (std::size_t ix = 0; ix < poll_results.size(); ++ix) {
               auto request = zmqSvc().receive<std::string>(*subscriber_ping_socket, zmq::recv_flags::dontwait);
               assert(request == "present");
               zmqSvc().send(*subscriber_ping_socket, "roger");
               ++N_subscribers_confirmed;
            }
         }
         zmqSvc().send(*mw_pub_, true);

      } else if (process_manager.is_queue()) {
         // first the queue-worker sockets
         // do resize instead of reserve so that the unique_ptrs are initialized
         // (to nullptr) so that we can do reset below, alternatively you can do
         // push/emplace_back with move or something
         qw_push_.resize(process_manager.N_workers());
         qw_pull_.resize(process_manager.N_workers());
         qw_push_poller_.resize(process_manager.N_workers());
         qw_pull_poller_.resize(process_manager.N_workers());
         for (std::size_t ix = 0; ix < process_manager.N_workers(); ++ix) {
            std::stringstream push_name, pull_name;
            // push
            qw_push_[ix].reset(zmqSvc().socket_ptr(zmq::socket_type::push));
            push_name << "ipc:///tmp/roofitMP_from_queue_to_worker_" << ix;
            qw_push_[ix]->bind(push_name.str());

            qw_push_poller_[ix].register_socket(*qw_push_[ix], zmq::event_flags::pollout);

            // pull
            qw_pull_[ix].reset(zmqSvc().socket_ptr(zmq::socket_type::pull));
            pull_name << "ipc:///tmp/roofitMP_from_worker_" << ix << "_to_queue";
            qw_pull_[ix]->bind(pull_name.str());

            qw_pull_poller_[ix].register_socket(*qw_pull_[ix], zmq::event_flags::pollin);
         }

         // then the master-queue sockets
         mq_push_.reset(zmqSvc().socket_ptr(zmq::socket_type::push));
         auto rc = zmq_setsockopt(*mq_push_, ZMQ_SNDHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         mq_push_->connect("ipc:///tmp/roofitMP_from_queue_to_master");

         mq_push_poller_.register_socket(*mq_push_, zmq::event_flags::pollout);

         mq_pull_.reset(zmqSvc().socket_ptr(zmq::socket_type::pull));
         rc = zmq_setsockopt(*mq_pull_, ZMQ_RCVHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         mq_pull_->connect("ipc:///tmp/roofitMP_from_master_to_queue");

         mq_pull_poller_.register_socket(*mq_pull_, zmq::event_flags::pollin);

         close_MQ_on_destruct_ = true;
         close_QW_container_on_destruct_ = true;
      } else if (process_manager.is_worker()) {
         // we only need one queue-worker pipe on the worker
         qw_push_poller_.resize(1);
         qw_pull_poller_.resize(1);

         std::stringstream push_name, pull_name;
         // push
         this_worker_qw_push_.reset(zmqSvc().socket_ptr(zmq::socket_type::push));
         push_name << "ipc:///tmp/roofitMP_from_worker_" << process_manager.worker_id() << "_to_queue";
         this_worker_qw_push_->connect(push_name.str());

         qw_push_poller_[0].register_socket(*this_worker_qw_push_, zmq::event_flags::pollout);

         // pull
         this_worker_qw_pull_.reset(zmqSvc().socket_ptr(zmq::socket_type::pull));
         pull_name << "ipc:///tmp/roofitMP_from_queue_to_worker_" << process_manager.worker_id();
         this_worker_qw_pull_->connect(pull_name.str());

         qw_pull_poller_[0].register_socket(*this_worker_qw_pull_, zmq::event_flags::pollin);

         mw_sub_.reset(zmqSvc().socket_ptr(zmq::socket_type::sub));
         auto rc = zmq_setsockopt(*mw_sub_, ZMQ_RCVHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         rc = zmq_setsockopt(*mw_sub_, ZMQ_SUBSCRIBE, "", 0);
         assert(rc == 0);
         mw_sub_->connect("ipc:///tmp/roofitMP_from_master_to_workers");
         mw_sub_poller_.register_socket(*mw_sub_, zmq::event_flags::pollin);

         wm_push_.reset(zmqSvc().socket_ptr(zmq::socket_type::push));
         rc = zmq_setsockopt(*wm_push_, ZMQ_SNDHWM, &hwm, sizeof hwm);
         assert(rc == 0);
         wm_push_->connect("ipc:///tmp/roofitMP_from_workers_to_master");

         // check publisher connection and then wait until all subscribers are connected
         ZmqLingeringSocketPtr<> subscriber_ping_socket{zmqSvc().socket_ptr(zmq::socket_type::req)};
         subscriber_ping_socket->connect("ipc:///tmp/roofitMP_subscriber_ping_socket");
         auto all_connected = zmqSvc().receive<bool>(*mw_sub_);
         zmqSvc().send(*subscriber_ping_socket, "present");
         auto reply = zmqSvc().receive<std::string>(*subscriber_ping_socket);
         assert(reply == "roger");

         while (!all_connected) {
            all_connected = zmqSvc().receive<bool>(*mw_sub_);
         }

         close_this_QW_on_destruct_ = true;
      } else {
         // should never get here
         throw std::runtime_error("Messenger ctor: I'm neither master, nor queue, nor a worker");
      }
   } catch (zmq::error_t &e) {
      std::cerr << e.what() << " -- errnum: " << e.num() << std::endl;
      throw;
   };
}

Messenger::~Messenger()
{
   printf("Messenger dtor on PID %d\n", getpid());
   if (close_MQ_on_destruct_) {
      try {
         mq_push_.reset(nullptr);
         mq_pull_.reset(nullptr);
         mw_pub_.reset(nullptr);
         wm_pull_.reset(nullptr);
      } catch (const std::exception &e) {
         std::cerr << "WARNING: something in Messenger dtor threw an exception! Original exception message:\n"
                   << e.what() << std::endl;
      }
   }
   if (close_this_QW_on_destruct_) {
      this_worker_qw_push_.reset(nullptr);
      this_worker_qw_pull_.reset(nullptr);
      mw_sub_.reset(nullptr);
      wm_push_.reset(nullptr);
   }
   if (close_QW_container_on_destruct_) {
      for (auto &socket : qw_push_) {
         socket.reset(nullptr);
      }
      for (auto &socket : qw_pull_) {
         socket.reset(nullptr);
      }
   }
   zmqSvc().close_context();
}

void Messenger::test_send(X2X ping_value, test_snd_pipes snd_pipe, std::size_t worker_id)
{
   try {
      switch (snd_pipe) {
      case test_snd_pipes::M2Q: {
         send_from_master_to_queue(ping_value);
         break;
      }
      case test_snd_pipes::Q2M: {
         send_from_queue_to_master(ping_value);
         break;
      }
      case test_snd_pipes::Q2W: {
         send_from_queue_to_worker(worker_id, ping_value);
         break;
      }
      case test_snd_pipes::W2Q: {
         send_from_worker_to_queue(ping_value);
         break;
      }
      }
   } catch (zmq::error_t &e) {
      if (e.num() == EAGAIN) {
         throw std::runtime_error("Messenger::test_connections: SEND over master-queue connection timed out!");
      } else {
         throw;
      }
   }
}

void Messenger::test_receive(X2X expected_ping_value, test_rcv_pipes rcv_pipe, std::size_t worker_id)
{
   X2X handshake = X2X::initial_value;

   std::size_t max_tries = 3, tries = 0;
   bool carry_on = true;
   while (carry_on && (tries++ < max_tries)) {
      try {
         switch (rcv_pipe) {
         case test_rcv_pipes::fromMonQ: {
            handshake = receive_from_master_on_queue<X2X>();
            break;
         }
         case test_rcv_pipes::fromQonM: {
            handshake = receive_from_queue_on_master<X2X>();
            break;
         }
         case test_rcv_pipes::fromQonW: {
            handshake = receive_from_queue_on_worker<X2X>();
            break;
         }
         case test_rcv_pipes::fromWonQ: {
            handshake = receive_from_worker_on_queue<X2X>(worker_id);
            break;
         }
         }
         carry_on = false;
      } catch (ZMQ::ppoll_error_t &e) {
         auto response = handle_zmq_ppoll_error(e);
         if (response == zmq_ppoll_error_response::abort) {
            throw std::runtime_error("EINTR in test_receive and SIGTERM received, aborting\n");
         } else if (response == zmq_ppoll_error_response::unknown_eintr) {
            printf("EINTR in test_receive but no SIGTERM received, try %zu\n", tries);
            continue;
         } else if (response == zmq_ppoll_error_response::retry) {
            printf("EAGAIN in test_receive, try %zu\n", tries);
            continue;
         }
      } catch (zmq::error_t &e) {
         if (e.num() == EAGAIN) {
            throw std::runtime_error("Messenger::test_connections: RECEIVE over master-queue connection timed out!");
         } else {
            printf("unhandled zmq::error_t (not a ppoll_error_t) in Messenger::test_receive with errno %d: %s\n",
                   e.num(), e.what());
            throw;
         }
      }
   }

   if (handshake != expected_ping_value) {
      throw std::runtime_error(
         "Messenger::test_connections: RECEIVE over master-queue connection failed, did not receive expected value!");
   }
}

/// \brief Test whether push-pull sockets are working
///
/// \note This function tests the PUSH-PULL socket pairs only. The PUB-SUB sockets are already tested in the
/// constructor.
///
/// \param process_manager ProcessManager object used to instantiate this object. Used to identify which process we are
/// running on and hence which sockets need to be tested.
void Messenger::test_connections(const ProcessManager &process_manager)
{
   // Before blocking SIGTERM, set the signal handler, so we can also check after blocking whether a signal occurred
   // In our case, we already set it in the ProcessManager after forking to the queue and worker processes.
   sigset_t sigmask;
   sigemptyset(&sigmask);
   sigaddset(&sigmask, SIGTERM);
   sigprocmask(SIG_BLOCK, &sigmask, &ppoll_sigmask);

   if (process_manager.is_master()) {
      test_send(X2X::ping, test_snd_pipes::M2Q, -1);
      test_receive(X2X::pong, test_rcv_pipes::fromQonM, -1);
      test_receive(X2X::ping, test_rcv_pipes::fromQonM, -1);
      test_send(X2X::pong, test_snd_pipes::M2Q, -1);
   } else if (process_manager.is_queue()) {
      ZeroMQPoller poller;
      std::size_t mq_index;
      std::tie(poller, mq_index) = create_queue_poller();

      for (std::size_t ix = 0; ix < process_manager.N_workers(); ++ix) {
         test_send(X2X::ping, test_snd_pipes::Q2W, ix);
      }

      while (!process_manager.sigterm_received() && (poller.size() > 0)) {
         // poll: wait until status change (-1: infinite timeout)
         std::vector<std::pair<size_t, zmq::event_flags>> poll_result;
         bool abort;
         std::tie(poll_result, abort) = careful_ppoll(poller, ppoll_sigmask);
         if (abort)
            break;

         // then process incoming messages from sockets
         for (auto readable_socket : poll_result) {
            // message comes from the master/queue socket (first element):
            if (readable_socket.first == mq_index) {
               test_receive(X2X::ping, test_rcv_pipes::fromMonQ, -1);
               test_send(X2X::pong, test_snd_pipes::Q2M, -1);
               test_send(X2X::ping, test_snd_pipes::Q2M, -1);
               test_receive(X2X::pong, test_rcv_pipes::fromMonQ, -1);
               printf("unregistering socket mq_pull on PID %d", getpid());
               poller.unregister_socket(*mq_pull_);
            } else { // from a worker socket
               // TODO: dangerous assumption for this_worker_id, may become invalid if we allow multiple queue_loops on
               // the same process!
               auto this_worker_id = readable_socket.first - 1; // TODO: replace with a more reliable lookup

               test_receive(X2X::pong, test_rcv_pipes::fromWonQ, this_worker_id);
               test_receive(X2X::ping, test_rcv_pipes::fromWonQ, this_worker_id);
               test_send(X2X::pong, test_snd_pipes::Q2W, this_worker_id);

               printf("unregistering socket qw_pull_[%lu] on PID %d", this_worker_id, getpid());
               poller.unregister_socket(*qw_pull_[this_worker_id]);
            }
         }
      }

   } else if (process_manager.is_worker()) {
      test_receive(X2X::ping, test_rcv_pipes::fromQonW, -1);
      test_send(X2X::pong, test_snd_pipes::W2Q, -1);
      test_send(X2X::ping, test_snd_pipes::W2Q, -1);
      test_receive(X2X::pong, test_rcv_pipes::fromQonW, -1);
   } else {
      // should never get here
      throw std::runtime_error("Messenger::test_connections: I'm neither master, nor queue, nor a worker");
   }

   // clean up signal management modifications
   sigprocmask(SIG_SETMASK, &ppoll_sigmask, nullptr);
   printf("done with test_connections on PID %d\n", getpid());
}

/// Helper function that creates a poller for Queue::loop()
std::pair<ZeroMQPoller, std::size_t> Messenger::create_queue_poller()
{
   ZeroMQPoller poller;
   std::size_t mq_index = poller.register_socket(*mq_pull_, zmq::event_flags::pollin);
   for (auto &s : qw_pull_) {
      poller.register_socket(*s, zmq::event_flags::pollin);
   }
   return {std::move(poller), mq_index};
}

/// Helper function that creates a poller for worker_loop()
std::pair<ZeroMQPoller, std::size_t> Messenger::create_worker_poller()
{
   ZeroMQPoller poller;
   poller.register_socket(*this_worker_qw_pull_, zmq::event_flags::pollin);
   std::size_t mw_sub_index = poller.register_socket(*mw_sub_, zmq::event_flags::pollin);
   return {std::move(poller), mw_sub_index};
}

// -- WORKER - QUEUE COMMUNICATION --

void Messenger::send_from_worker_to_queue() {}

void Messenger::send_from_queue_to_worker(std::size_t /*this_worker_id*/) {}

// -- QUEUE - MASTER COMMUNICATION --

void Messenger::send_from_queue_to_master() {}

void Messenger::send_from_master_to_queue() {}

/// Set the flag used in all send functions; 0, ZMQ_DONTWAIT, ZMQ_SNDMORE or bitwise combination
void Messenger::set_send_flag(zmq::send_flags flag)
{
   send_flag_ = flag;
}

// -- MASTER - WORKER COMMUNICATION --

void Messenger::publish_from_master_to_workers() {}

void Messenger::send_from_worker_to_master() {}

// for debugging
#define PROCESS_VAL(p) \
   case (p): s = #p; break;

std::ostream &operator<<(std::ostream &out, const M2Q value)
{
   std::string s;
   switch (value) {
      PROCESS_VAL(M2Q::enqueue);
   default: s = std::to_string(static_cast<int>(value));
   }
   return out << s;
}

std::ostream &operator<<(std::ostream &out, const W2Q value)
{
   std::string s;
   switch (value) {
      PROCESS_VAL(W2Q::dequeue);
   default: s = std::to_string(static_cast<int>(value));
   }
   return out << s;
}

std::ostream &operator<<(std::ostream &out, const Q2W value)
{
   std::string s;
   switch (value) {
      PROCESS_VAL(Q2W::dequeue_rejected);
      PROCESS_VAL(Q2W::dequeue_accepted);
   default: s = std::to_string(static_cast<int>(value));
   }
   return out << s;
}

std::ostream &operator<<(std::ostream &out, const X2X value)
{
   std::string s;
   switch (value) {
      PROCESS_VAL(X2X::ping);
      PROCESS_VAL(X2X::pong);
   default: s = std::to_string(static_cast<int>(value));
   }
   return out << s;
}

#undef PROCESS_VAL

/// Function called from send and receive template functions in debug builds
/// used to monitor the messages that are going to be sent or are received.
/// By defining this in the implementation file, compilation is a lot faster
/// during debugging of Messenger or communication protocols.
void Messenger::debug_print(std::string s)
{
   // print 's' when debugging
   std::cerr << s << '\n';
}

} // namespace MultiProcess
} // namespace RooFit
