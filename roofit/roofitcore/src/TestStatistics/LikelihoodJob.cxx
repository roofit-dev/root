/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include "LikelihoodJob.h"

#include "LikelihoodSerial.h"
#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/ProcessManager.h"
#include "RooFit/MultiProcess/Queue.h"
#include "RooFit/MultiProcess/Job.h"
#include "RooFit/MultiProcess/types.h"
#include "RooFit/MultiProcess/Config.h"
#include "RooFit/MultiProcess/ProcessTimer.h"
#include "RooFit/TestStatistics/RooAbsL.h"
#include "RooFit/TestStatistics/RooUnbinnedL.h"
#include "RooFit/TestStatistics/RooBinnedL.h"
#include "RooFit/TestStatistics/RooSubsidiaryL.h"
#include "RooFit/TestStatistics/RooSumL.h"
#include "RooRealVar.h"
#include "TClass.h"
#include <ctime>
#include <chrono>
#include <iostream>

#include <sys/resource.h>
#include <sys/time.h>

namespace RooFit {
namespace TestStatistics {

LikelihoodJob::LikelihoodJob(std::shared_ptr<RooAbsL> likelihood,
                             std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean)
   : LikelihoodJob(std::move(likelihood), std::move(calculation_is_clean),
                   std::make_shared<std::vector<ROOT::Math::KahanSum<double>>>(), std::make_shared<std::vector<ROOT::Math::KahanSum<double>>>())
{
}

LikelihoodJob::LikelihoodJob(std::shared_ptr<RooAbsL> likelihood,
                             std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean,
                             std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> offsets,
                             std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> offsets_save)
   : LikelihoodWrapper(std::move(likelihood), std::move(calculation_is_clean), std::move(offsets),
                       std::move(offsets_save)),
     n_event_tasks_(MultiProcess::Config::LikelihoodJob::defaultNEventTasks),
     n_component_tasks_(MultiProcess::Config::LikelihoodJob::defaultNComponentTasks),
     likelihood_serial_(std::make_unique<LikelihoodSerial>(likelihood_, calculation_is_clean_, component_offsets_, component_offsets_save_))
{
   init_vars();
   offsets_previous_ = *component_offsets_;
}

LikelihoodJob::LikelihoodJob(const LikelihoodJob &other)
   : MultiProcess::Job(other), LikelihoodWrapper(other), result_(other.result_), results_(other.results_), vars_(other.vars_),
     save_vars_(other.save_vars_), n_tasks_at_workers_(other.n_tasks_at_workers_), n_event_tasks_(other.n_event_tasks_),
     n_component_tasks_(other.n_component_tasks_), offsets_previous_(other.offsets_previous_),
     likelihood_serial_(other.likelihood_serial_->clone())
{
}

LikelihoodJob *LikelihoodJob::clone() const
{
   return new LikelihoodJob(*this);
}

// This is a separate function (instead of just in ctor) for historical reasons.
// Its predecessor RooRealMPFE::initVars() was used from multiple ctors, but also
// from RooRealMPFE::constOptimizeTestStatistic at the end, which makes sense,
// because it might change the set of variables. We may at some point want to do
// this here as well.
void LikelihoodJob::init_vars()
{
   // Empty current lists
   vars_.removeAll();
   save_vars_.removeAll();

   // Retrieve non-constant parameters
   auto vars = std::make_unique<RooArgSet>(
      *likelihood_->getParameters()); // TODO: make sure this is the right list of parameters, compare to original
                                      // implementation in RooRealMPFE.cxx
   RooArgList varList(*vars);

   // Save in lists
   vars_.add(varList);
   save_vars_.addClone(varList);

   if (MultiProcess::Config::isInLinesearch_)
   {
//    std::cout << "pruning client list" << std::endl;
      for (auto& var: vars_) var->pruneClientValueListHack();
   }
}

void LikelihoodJob::update_state()
{
   if (get_manager()->process_manager().is_worker()) {
//      std::chrono::time_point<std::chrono::steady_clock> worker_update_state_begin = std::chrono::steady_clock::now();
      bool more;

//      std::chrono::time_point<std::chrono::steady_clock> worker_update_state_1 = std::chrono::steady_clock::now();
      auto mode = get_manager()->messenger().receive_from_master_on_worker<update_state_mode>(&more);
      assert(more);
//      std::cout << "in linesearch calculation " <<  MultiProcess::Config::isInLinesearch_ <<std::endl;
      if (MultiProcess::Config::isInLinesearch_ & !MultiProcess::Config::clientListPruned_)
      {
//         std::cout << "pruning client list" << std::endl;
         for (auto& var: vars_) var->pruneClientValueListHack();

         MultiProcess::Config::clientListPruned_ = true;
      }
 //     if (MultiProcess::Config::isInLinesearch_) std::cout << "worker update state 1: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_1).count() << std::endl;


      switch (mode) {
      case update_state_mode::parameters: {
//      std::chrono::time_point<std::chrono::steady_clock> worker_update_state_2 = std::chrono::steady_clock::now();
         state_id_ = get_manager()->messenger().receive_from_master_on_worker<RooFit::MultiProcess::State>(&more);
         assert(more);
//      if (MultiProcess::Config::isInLinesearch_) std::cout << "worker update state 2: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_2).count() << std::endl;
//      std::chrono::time_point<std::chrono::steady_clock> worker_update_state_3 = std::chrono::steady_clock::now();
         auto message = get_manager()->messenger().receive_from_master_on_worker<zmq::message_t>(&more);
//      if (MultiProcess::Config::isInLinesearch_) std::cout << "worker update state 3: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_3).count() << std::endl;
//      std::chrono::time_point<std::chrono::steady_clock> worker_update_state_4 = std::chrono::steady_clock::now();
         auto message_begin = message.data<update_state_t>();
         auto message_end = message_begin + message.size() / sizeof(update_state_t);
         std::vector<update_state_t> to_update(message_begin, message_end);
//      if (MultiProcess::Config::isInLinesearch_) std::cout << "worker update state 4: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_4).count() << std::endl;
//      std::chrono::time_point<std::chrono::steady_clock> worker_update_state_5 = std::chrono::steady_clock::now();
      // if (MultiProcess::Config::isInLinesearch_) std::cout << "starting loop!!" << std::endl;
         for (auto const &item : to_update) {
            // std::chrono::time_point<std::chrono::steady_clock> loop_1 = std::chrono::steady_clock::now();
            RooRealVar *rvar = (RooRealVar *)vars_.at(item.var_index);
            // if (MultiProcess::Config::isInLinesearch_) std::cout << "loop 1: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - loop_1).count() << std::endl;
            // std::chrono::time_point<std::chrono::steady_clock> loop_2 = std::chrono::steady_clock::now();
            rvar->setVal(static_cast<double>(item.value));
            //for (auto& client: rvar->valueClients())
            //   std::cout << client->IsA()->GetName() << "::"<< client->GetName() << "::" << client->operMode() << std::endl;
            // if (MultiProcess::Config::isInLinesearch_) std::cout << "loop 2: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - loop_2).count() << std::endl;
            // std::chrono::time_point<std::chrono::steady_clock> loop_3 = std::chrono::steady_clock::now();
            if (rvar->isConstant() != item.is_constant) {
               rvar->setConstant(static_cast<bool>(item.is_constant));
            }
            // if (MultiProcess::Config::isInLinesearch_) std::cout << "loop 3: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - loop_3).count() << std::endl;
         }
      // if (MultiProcess::Config::isInLinesearch_) std::cout << "ending loop" << std::endl;
//      if (MultiProcess::Config::isInLinesearch_) std::cout << "worker update state 5: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_5).count() << std::endl;
//      std::chrono::time_point<std::chrono::steady_clock> worker_update_state_6 = std::chrono::steady_clock::now();
         if (more) {
            // offsets also incoming
            auto offsets_message = get_manager()->messenger().receive_from_master_on_worker<zmq::message_t>(&more);
            assert(!more);
            auto offsets_message_begin = offsets_message.data<ROOT::Math::KahanSum<double>>();
            std::size_t N_offsets = offsets_message.size() / sizeof(ROOT::Math::KahanSum<double>);
            component_offsets_->reserve(N_offsets);
            auto offsets_message_end = offsets_message_begin + N_offsets;
            std::copy(offsets_message_begin, offsets_message_end, component_offsets_->begin());
         }
//      if (MultiProcess::Config::isInLinesearch_) std::cout << "worker update state 6: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_6).count() << std::endl;
         break;
      }
      case update_state_mode::offsetting: {
         LikelihoodWrapper::enableOffsetting(get_manager()->messenger().receive_from_master_on_worker<bool>(&more));
         assert(!more);
         break;
      }
      }
//      if (MultiProcess::Config::isInLinesearch_) std::cout << "worker update state: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - worker_update_state_begin).count() << std::endl;
   }
}

/// \warning In automatic mode, this function can start MultiProcess (forks, starts workers, etc)!
std::size_t LikelihoodJob::getNEventTasks()
{
   std::size_t val = n_event_tasks_;
   if (val == MultiProcess::Config::LikelihoodJob::automaticNEventTasks) {
      val = get_manager()->process_manager().N_workers();
   }
   if (val > likelihood_->getNEvents()) {
      val = likelihood_->getNEvents();
   }
   return val;
}


std::size_t LikelihoodJob::getNComponentTasks()
{
   std::size_t val = n_component_tasks_;
   if (val == MultiProcess::Config::LikelihoodJob::automaticNComponentTasks) {
      val = 1;
   }
   if (val > likelihood_->getNComponents()) {
      val = likelihood_->getNComponents();
   }
   return val;
}

void LikelihoodJob::updateWorkersParameters()
{
   if (get_manager()->process_manager().is_master()) {
      bool valChanged = false;
      bool constChanged = false;
      std::vector<update_state_t> to_update;
      for (std::size_t ix = 0u; ix < static_cast<std::size_t>(vars_.getSize()); ++ix) {
         valChanged = !vars_[ix].isIdentical(save_vars_[ix], true);
         constChanged = (vars_[ix].isConstant() != save_vars_[ix].isConstant());

         if (valChanged || constChanged) {
            if (constChanged) {
               ((RooRealVar *)&save_vars_[ix])->setConstant(vars_[ix].isConstant());
            }
            // TODO: Check with Wouter why he uses copyCache in MPFE; makes it very difficult to extend, because
            // copyCache is protected (so must be friend). Moved setting value to if-block below.
            //          _saveVars[ix].copyCache(&_vars[ix]);

            // send message to queue (which will relay to workers)
            RooAbsReal *rar_val = dynamic_cast<RooAbsReal *>(&vars_[ix]);
            if (rar_val) {
               double val = rar_val->getVal();
               dynamic_cast<RooRealVar *>(&save_vars_[ix])->setVal(val);
               bool isC = vars_[ix].isConstant();
               to_update.push_back(update_state_t{ix, val, isC});
            }
         }
      }
      bool update_offsets = isOffsetting() && component_offsets_ != nullptr && *component_offsets_ != offsets_previous_;
      if (!to_update.empty() || update_offsets) {
         ++state_id_;
         zmq::message_t message(to_update.begin(), to_update.end());
         // always send Job id first! This is used in worker_loop to route the
         // update_state call to the correct Job.
         if (update_offsets) {
            zmq::message_t offsets_message(component_offsets_->begin(), component_offsets_->end());
            get_manager()->messenger().publish_from_master_to_workers(id_, update_state_mode::parameters, state_id_, std::move(message), std::move(offsets_message));
            offsets_previous_ = *component_offsets_;
         } else {
            get_manager()->messenger().publish_from_master_to_workers(id_, update_state_mode::parameters, state_id_, std::move(message));
         }
      }
   }
}

void LikelihoodJob::updateWorkersOffsetting()
{
   get_manager()->messenger().publish_from_master_to_workers(id_, update_state_mode::offsetting, isOffsetting());
}

void LikelihoodJob::evaluate()
{
   if (get_manager()->process_manager().is_master()) {
      // evaluate the serial likelihood to set the offsets
//      std::chrono::time_point<std::chrono::steady_clock> offsets_begin = std::chrono::steady_clock::now();
      if (do_offset_ && component_offsets_->empty()) {
         likelihood_serial_->evaluate();
         std::cout << "evaluating serial likelihood on master" << std::endl;
         // note: we don't need to get the offsets from the serial likelihood, because they are already coupled through
         // the shared_ptr
      }
//      if (MultiProcess::Config::isInLinesearch_) std::cout << "master offsets: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - offsets_begin).count() << std::endl;

//      std::chrono::time_point<std::chrono::steady_clock> update_workers_parameters_begin = std::chrono::steady_clock::now();

      // update parameters that changed since last calculation (or creation if first time)
      updateWorkersParameters();

//      if (MultiProcess::Config::isInLinesearch_) std::cout << "master update workers parameters: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - update_workers_parameters_begin).count() << std::endl;

//      std::chrono::time_point<std::chrono::steady_clock> fill_receive_tasks_begin = std::chrono::steady_clock::now();

      // master fills queue with tasks
      auto N_tasks = getNEventTasks() * getNComponentTasks();
      for (std::size_t ix = 0; ix < N_tasks; ++ix) {
         get_manager()->queue()->add({id_, state_id_, ix});
      }
      n_tasks_at_workers_ = N_tasks;

      // wait for task results back from workers to master
      gather_worker_results();

//      if (MultiProcess::Config::isInLinesearch_) std::cout << "master fill receive tasks: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - fill_receive_tasks_begin).count() << std::endl;

//      std::chrono::time_point<std::chrono::steady_clock> process_results_begin = std::chrono::steady_clock::now();

      // Note: initializing result_ to results_[0] instead of zero-initializing it makes
      // a difference due to Kahan sum precision. This way, a single-worker run gives
      // the same result as a run with serial likelihood. Adding the terms to a zero
      // initial sum can cancel the carry in some cases, causing divergent values.
      result_ = results_[0];
      for (auto item_it = results_.cbegin() + 1; item_it != results_.cend(); ++item_it) {
         result_ += *item_it;
      }
      results_.clear();
      // if (MultiProcess::Config::isInLinesearch_) std::cout << "master process results: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - process_results_begin).count() << std::endl;
   }
}

// --- RESULT LOGISTICS ---

void LikelihoodJob::send_back_task_result_from_worker(std::size_t /*task*/)
{
   task_result_t task_result{id_, result_.Result(), result_.Carry()};
   zmq::message_t message(sizeof(task_result_t));
   memcpy(message.data(), &task_result, sizeof(task_result_t));
   get_manager()->messenger().send_from_worker_to_master(std::move(message));
}

bool LikelihoodJob::receive_task_result_on_master(const zmq::message_t &message)
{
   auto task_result = message.data<task_result_t>();
   results_.emplace_back(task_result->value, task_result->carry);
   --n_tasks_at_workers_;
   bool job_completed = (n_tasks_at_workers_ == 0);
   return job_completed;
}

// --- END OF RESULT LOGISTICS ---

void LikelihoodJob::evaluate_task(std::size_t task)
{
//   struct rusage usage;
//  struct timeval start, end;
//   std::chrono::time_point<std::chrono::steady_clock> eval_task_begin = std::chrono::steady_clock::now();

//   getrusage(RUSAGE_SELF, &usage);
//   start = usage.ru_utime;

   assert(get_manager()->process_manager().is_worker());

   double section_first = 0;
   double section_last = 1;
   if (getNEventTasks() > 1) {
      std::size_t event_task = task % getNEventTasks();
      std::size_t N_events = likelihood_->numDataEntries();
      if (event_task > 0) {
         std::size_t first = N_events * event_task / getNEventTasks();
         section_first = static_cast<double>(first) / N_events;
      }
      if (event_task < getNEventTasks() - 1) {
         std::size_t last = N_events * (event_task + 1) / getNEventTasks();
         section_last = static_cast<double>(last) / N_events;
      }
   }

   switch (likelihood_type_) {
   case LikelihoodType::unbinned:
   case LikelihoodType::binned: {
      result_ = likelihood_->evaluatePartition({section_first, section_last}, 0, 0);
      if (do_offset_ && section_last == 1) {
         // we only subtract at the end of event sections, otherwise the offset is subtracted for each event split
         result_ -= (*component_offsets_)[0];
      }
      break;
   }
   case LikelihoodType::subsidiary: {
      result_ = likelihood_->evaluatePartition({0, 1}, 0, 0);
      if (do_offset_ && offsetting_mode_ == OffsettingMode::full) {
         result_ -= (*component_offsets_)[0];
      }
      break;
   }
   case LikelihoodType::sum: {
      std::size_t components_first = 0;
      std::size_t components_last = likelihood_->getNComponents();
      if (getNComponentTasks() > 1) {
         std::size_t component_task = task / getNEventTasks();
         components_first = likelihood_->getNComponents() * component_task / getNComponentTasks();
         if (component_task == getNComponentTasks() - 1) {
            components_last = likelihood_->getNComponents();
         } else {
            components_last = likelihood_->getNComponents() * (component_task + 1) / getNComponentTasks();
         }
      }

      result_ = ROOT::Math::KahanSum<double>();
      for (std::size_t comp_ix = components_first; comp_ix < components_last; ++comp_ix) {
         auto component_result = likelihood_->evaluatePartition({section_first, section_last}, comp_ix, comp_ix + 1);
         if (do_offset_ && section_last == 1 && (*component_offsets_)[comp_ix] != ROOT::Math::KahanSum<double>(0, 0)) {
            // we only subtract at the end of event sections, otherwise the offset is subtracted for each event split
            result_ += (component_result - (*component_offsets_)[comp_ix]);
         } else {
            result_ += component_result;
         }
      }

      break;
   }
   }

//  getrusage(RUSAGE_SELF, &usage);
//  end = usage.ru_utime;

//   if (MultiProcess::Config::isInLinesearch_) std::cout << "evaluate task (cpu): " << end.tv_usec - start.tv_usec << " | ";
//   if (MultiProcess::Config::isInLinesearch_) std::cout << "(wtime): " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - eval_task_begin).count() << std::endl;
}

void LikelihoodJob::enableOffsetting(bool flag)
{
   likelihood_serial_->enableOffsetting(flag);
   LikelihoodWrapper::enableOffsetting(flag);
   if (RooFit::MultiProcess::JobManager::is_instantiated()) {
      printf("WARNING: when calling MinuitFcnGrad::setOffsetting after the run has already been started the MinuitFcnGrad::likelihood_in_gradient object (a LikelihoodSerial) on the workers can no longer be updated! This function (LikelihoodJob::enableOffsetting) can in principle be used outside of MinuitFcnGrad, but be aware of this limitation. To do a minimization with a different offsetting setting, please delete all RooFit::MultiProcess based objects so that the forked processes are killed and then set up a new RooMinimizer.\n");
      updateWorkersOffsetting();
   }
}

#define PROCESS_VAL(p) \
   case (p): s = #p; break;

std::ostream &operator<<(std::ostream &out, const LikelihoodJob::update_state_mode value)
{
   std::string s;
   switch (value) {
      PROCESS_VAL(LikelihoodJob::update_state_mode::offsetting);
      PROCESS_VAL(LikelihoodJob::update_state_mode::parameters);
   default: s = std::to_string(static_cast<int>(value));
   }
   return out << s;
}

#undef PROCESS_VAL

} // namespace TestStatistics
} // namespace RooFit
