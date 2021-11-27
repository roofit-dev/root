/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2020, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#include "TestStatistics/LikelihoodGradientJob.h"

#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/Messenger.h"
#include "RooFit/MultiProcess/Queue.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"

#include "Minuit2/MnStrategy.h"

namespace RooFit {
namespace TestStatistics {

LikelihoodGradientJob::LikelihoodGradientJob(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean, std::size_t N_dim,
                                             RooMinimizer *minimizer)
   : LikelihoodGradientWrapper(std::move(likelihood), std::move(calculation_is_clean), N_dim, minimizer), grad_(N_dim)
{
   // Note to future maintainers: take care when storing the minimizer_fcn pointer. The
   // RooAbsMinimizerFcn subclasses may get cloned inside MINUIT, which means the pointer
   // should also somehow be updated in this class.
//   N_tasks = minimizer_fcn->get_nDim();
   N_tasks_ = N_dim;
   minuit_internal_x_.reserve(N_dim);
   // TODO: make sure that the full gradients are sent back so that the
   // derivator will depart from correct state next step everywhere!
}

LikelihoodGradientJob::LikelihoodGradientJob(const LikelihoodGradientJob &other)
   : LikelihoodGradientWrapper(other), grad_(other.grad_), gradf_(other.gradf_), N_tasks_(other.N_tasks_),
     minuit_internal_x_(other.minuit_internal_x_)
{
}

LikelihoodGradientJob *LikelihoodGradientJob::clone() const
{
   return new LikelihoodGradientJob(*this);
}

void LikelihoodGradientJob::synchronizeParameterSettings(
   ROOT::Math::IMultiGenFunction *function, const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings)
{
   gradf_.SetInitialGradient(function, parameter_settings, grad_);
}

void LikelihoodGradientJob::synchronizeWithMinimizer(const ROOT::Math::MinimizerOptions &options)
{
   setStrategy(options.Strategy());
   setErrorLevel(options.ErrorDef());
}

void LikelihoodGradientJob::setStrategy(int istrat)
{
   assert(istrat >= 0);
   ROOT::Minuit2::MnStrategy strategy(static_cast<unsigned int>(istrat));

   setStepTolerance(strategy.GradientStepTolerance());
   setGradTolerance(strategy.GradientTolerance());
   setNCycles(strategy.GradientNCycles());
}

void LikelihoodGradientJob::setStepTolerance(double step_tolerance) const
{
   gradf_.SetStepTolerance(step_tolerance);
}

void LikelihoodGradientJob::setGradTolerance(double grad_tolerance) const
{
   gradf_.SetGradTolerance(grad_tolerance);
}

void LikelihoodGradientJob::setNCycles(unsigned int ncycles) const
{
   gradf_.SetNCycles(ncycles);
}

void LikelihoodGradientJob::setErrorLevel(double error_level) const
{
   gradf_.SetErrorLevel(error_level);
}

///////////////////////////////////////////////////////////////////////////////
/// Job overrides:

void LikelihoodGradientJob::evaluate_task(std::size_t task)
{
   run_derivator(task);
}

// SYNCHRONIZATION FROM WORKERS TO MASTER

bool LikelihoodGradientJob::receive_task_result_on_master(const zmq::message_t & message)
{
   auto result = message.data<task_result_t>();
   grad_[result->task_id] = result->grad;
   --N_tasks_at_workers_;
   bool job_completed = (N_tasks_at_workers_ == 0);
   return job_completed;
}


// END SYNCHRONIZATION FROM WORKERS TO MASTER

///////////////////////////////////////////////////////////////////////////////
/// Calculation stuff (mostly duplicates of RooGradMinimizerFcn code):

void LikelihoodGradientJob::run_derivator(unsigned int i_component) const
{
   // Calculate the derivative etc for these parameters
//   auto parameter_values = _minimizer->get_function_parameter_values();
   grad_[i_component] = gradf_.FastPartialDerivative(
      minimizer_->getMultiGenFcn(), minimizer_->fitter()->Config().ParamsSettings(), i_component, grad_[i_component]);
}

///////////////////////////////////////////////////////////////////////////////
/// copy pasted and adapted from old MP::GradMinimizerFcn:

//void LikelihoodGradientJob::update_workers_state()
//{
//   auto get_time = []() {
//      return std::chrono::duration_cast<std::chrono::nanoseconds>(
//         std::chrono::high_resolution_clock::now().time_since_epoch())
//         .count();
//   };
//   decltype(get_time()) t1, t2;
//   t1 = get_time();
//   // TODO optimization: only send changed parameters (now sending all)
//   std::size_t ix = 0;
//   RooFit::MultiProcess::M2Q msg = RooFit::MultiProcess::M2Q::update_real;
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix, _grad[ix].derivative, false);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix + 1 * _minimizer->getNPar(), _grad[ix].second_derivative,
//                                                           false);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix + 2 * _minimizer->getNPar(), _grad[ix].step_size,
//                                                           false);
//   }
//
////   ix = 0;
////   auto parameter_values = _minimizer->get_function_parameter_values();
////   for (auto &parameter : parameter_values) {
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix + 3 * _minimizer->getNPar(), /*parameter*/ minuit_internal_x_[ix], false);
////      ++ix;
//   }
//   t2 = get_time();
//   printf("timestamps LikelihoodGradientJob::update_workers_state: %lld %lld\n", t1, t2);
//}

//void LikelihoodGradientJob::update_workers_state()
//{
//   // TODO optimization: only send changed parameters (now sending all)
//   std::size_t ix;
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(_grad[ix].derivative);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(_grad[ix].second_derivative);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(_grad[ix].step_size);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(minuit_internal_x_[ix]);
//   }
//}

void LikelihoodGradientJob::update_workers_state()
{
   // TODO optimization: only send changed parameters (now sending all)
   get_manager()->messenger().publish_from_master_to_workers(id_);
   zmq::message_t gradient_message(grad_.begin(), grad_.end());
   get_manager()->messenger().publish_from_master_to_workers(std::move(gradient_message));
   zmq::message_t minuit_internal_x_message(minuit_internal_x_.begin(), minuit_internal_x_.end());
   get_manager()->messenger().publish_from_master_to_workers(std::move(minuit_internal_x_message));
}


void LikelihoodGradientJob::calculate_all()
{
//   auto get_time = []() {
//      return std::chrono::duration_cast<std::chrono::nanoseconds>(
//                std::chrono::high_resolution_clock::now().time_since_epoch())
//         .count();
//   };
//   decltype(get_time()) t1, t2;

//   std::cout << "BABBELBOX" << std::endl;

   if (get_manager()->process_manager().is_master()) {
//      std::cout << "HAAAAAA" << std::endl;
      // do Grad, G2 and Gstep here and then just return results from the
      // separate functions below

      // update parameters and object states that changed since last calculation (or creation if first time)
//      t1 = get_time();
      update_workers_state();
//      t2 = get_time();

//      printf("wallclock [master] update_workers_state: %f\n", (t2 - t1) / 1.e9);

//      t1 = get_time();
      // master fills queue with tasks
      for (std::size_t ix = 0; ix < N_tasks_; ++ix) {
         MultiProcess::JobTask job_task(id_, ix);
         get_manager()->queue().add(job_task);
      }
      N_tasks_at_workers_ = N_tasks_;
//      t2 = get_time();

//      printf("wallclock [master] put job tasks in queue: %f\n", (t2 - t1) / 1.e9);

      // wait for task results back from workers to master (put into _grad)
//      t1 = get_time();
      gather_worker_results();
//      t2 = get_time();

//      printf("wallclock [master] gather_worker_results: %f\n", (t2 - t1) / 1.e9);

      calculation_is_clean_->gradient = true;
//      calculation_is_clean_->g2 = true;
//      calculation_is_clean_->gstep = true;
   }
}

void LikelihoodGradientJob::fillGradient(double *grad)
{
   if (get_manager()->process_manager().is_master()) {
      if (!calculation_is_clean_->gradient) {
         calculate_all();
      }

      // TODO: maybe make a flag to avoid this copy operation, but maybe not worth the effort
      // put the results from _grad into *grad
      for (Int_t ix = 0; ix < minimizer_->getNPar(); ++ix) {
         grad[ix] = grad_[ix].derivative;
      }
   }
}

void LikelihoodGradientJob::updateMinuitInternalParameterValues(const std::vector<double>& minuit_internal_x)
{
   minuit_internal_x_ = minuit_internal_x;
}

bool LikelihoodGradientJob::usesMinuitInternalValues()
{
   return true;
}

void LikelihoodGradientJob::update_state()
{
//   std::size_t ix;
//
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      _grad[ix].derivative = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      _grad[ix].second_derivative = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      _grad[ix].step_size = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      minuit_internal_x_[ix] = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }

   auto gradient_message = get_manager()->messenger().receive_from_master_on_worker<zmq::message_t>();
   auto gradient_message_begin = gradient_message.data<ROOT::Minuit2::DerivatorElement>();
   auto gradient_message_end = gradient_message_begin + gradient_message.size()/sizeof(ROOT::Minuit2::DerivatorElement);
   std::copy(gradient_message_begin, gradient_message_end, grad_.begin());

   auto minuit_internal_x_message = get_manager()->messenger().receive_from_master_on_worker<zmq::message_t>();
   auto minuit_internal_x_message_begin = minuit_internal_x_message.data<double>();
   auto minuit_internal_x_message_end = minuit_internal_x_message_begin + minuit_internal_x_message.size()/sizeof(double);
   std::copy(minuit_internal_x_message_begin, minuit_internal_x_message_end, minuit_internal_x_.begin());

   gradf_.SetupDifferentiate(minimizer_->getMultiGenFcn(), minuit_internal_x_.data(),
                             minimizer_->fitter()->Config().ParamsSettings());
}


} // namespace TestStatistics
} // namespace RooFit
