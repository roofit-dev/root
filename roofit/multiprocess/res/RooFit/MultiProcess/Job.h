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
#ifndef ROOT_ROOFIT_MultiProcess_Job_decl
#define ROOT_ROOFIT_MultiProcess_Job_decl

#include <string>
#include <zmq.hpp>

namespace RooFit {
namespace MultiProcess {

// forward declaration
class JobManager;

class Job {
public:
   explicit Job();
   Job(const Job &other);

   ~Job();

   virtual void evaluate_task(std::size_t task) = 0;
   virtual void update_state();

   virtual void send_back_task_result_from_worker(std::size_t task) = 0;
   virtual bool receive_task_result_on_master(const zmq::message_t &message) = 0;

   void gather_worker_results();

protected:
   JobManager *get_manager();

   std::size_t id_;

private:
   // do not use _manager directly, it must first be initialized! use get_manager()
   JobManager *_manager = nullptr;
};

} // namespace MultiProcess
} // namespace RooFit

#endif // ROOT_ROOFIT_MultiProcess_Job_decl
