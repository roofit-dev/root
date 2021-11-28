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
#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodJob
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodJob

#include "RooFit/MultiProcess/Job.h"
#include "RooFit/MultiProcess/types.h"
#include "TestStatistics/LikelihoodWrapper.h"
#include "RooArgList.h"

#include "Math/MinimizerOptions.h"

#include <vector>

namespace RooFit {
namespace TestStatistics {

class LikelihoodJob : public MultiProcess::Job, public LikelihoodWrapper {
public:
   LikelihoodJob(std::shared_ptr<RooAbsL> _likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean/*, RooMinimizer *minimizer*/);
   LikelihoodJob* clone() const override;

   void init_vars();

   void evaluate() override;
   inline ROOT::Math::KahanSum<double> getResult() const override { return result; }

   void updateWorkersParameters();  // helper for evaluate
   void updateWorkersOffsetting();  // helper for enableOffsetting

   // Job overrides:
   void evaluate_task(std::size_t task) override;
   void update_state() override;

   struct update_state_t {
      std::size_t var_index;
      double value;
      bool is_constant;
   };
   enum class update_state_mode : int {parameters, offsetting};

   // --- RESULT LOGISTICS ---
   struct task_result_t {
      std::size_t job_id; // job ID must always be the first part of any result message/type
      double value;
      double carry;
   };

   void send_back_task_result_from_worker(std::size_t task) override;
   bool receive_task_result_on_master(const zmq::message_t & message) override;

   void enableOffsetting(bool flag) override;

private:
   ROOT::Math::KahanSum<double> result;
   std::vector<ROOT::Math::KahanSum<double>> results;

   RooArgList vars_;      // Variables
   RooArgList save_vars_;  // Copy of variables

   LikelihoodType likelihood_type_;
   std::size_t N_tasks_at_workers_ = 0;
};

std::ostream &operator<<(std::ostream &out, const LikelihoodJob::update_state_mode value);

}
}

#endif // ROOT_ROOFIT_LikelihoodJob
