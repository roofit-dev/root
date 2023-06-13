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

#include <utility>
#include "RooFit/MultiProcess/ProcessManager.h"
#include "RooFit/MultiProcess/Queue.h"
#include "RooFit/MultiProcess/Job.h"
#include "RooFit/MultiProcess/types.h"
#include "RooFit/MultiProcess/Config.h"
#include "RooFit/TestStatistics/RooAbsL.h"
#include "RooFit/TestStatistics/RooUnbinnedL.h"
#include "RooFit/TestStatistics/RooBinnedL.h"
#include "RooFit/TestStatistics/RooSubsidiaryL.h"
#include "RooFit/TestStatistics/RooSumL.h"
#include "RooRealVar.h"
#include "RooNaNPacker.h"

#include "TMath.h" // IsNaN

namespace RooFit {
namespace TestStatistics {

LikelihoodJob::LikelihoodJob(std::shared_ptr<RooAbsL> likelihood,
                             std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean)
   : LikelihoodJob(std::move(likelihood), std::move(calculation_is_clean),
                   std::make_shared<std::vector<ROOT::Math::KahanSum<double>>>(),
                   std::make_shared<std::vector<ROOT::Math::KahanSum<double>>>())
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
     likelihood_serial_(std::make_unique<LikelihoodSerial>(likelihood_, calculation_is_clean_, component_offsets_,
                                                           component_offsets_save_))
{
   init_vars();
   offsets_previous_ = *component_offsets_;
}

LikelihoodJob::LikelihoodJob(const LikelihoodJob &other)
   : MultiProcess::Job(other), LikelihoodWrapper(other), result_(other.result_), results_(other.results_),
     vars_(other.vars_), save_vars_(other.save_vars_), n_tasks_at_workers_(other.n_tasks_at_workers_),
     n_event_tasks_(other.n_event_tasks_), n_component_tasks_(other.n_component_tasks_),
     offsets_previous_(other.offsets_previous_), likelihood_serial_(other.likelihood_serial_->clone())
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
   std::unique_ptr<RooArgSet> vars{likelihood_->getParameters()};
   // TODO: make sure this is the right list of parameters, compare to original
   // implementation in RooRealMPFE.cxx

   RooArgList varList(*vars);

   // Save in lists
   vars_.add(varList);
   save_vars_.addClone(varList);
}

void LikelihoodJob::update_state()
{
   if (get_manager()->process_manager().is_worker()) {
      bool more;

      auto mode = get_manager()->messenger().receive_from_master_on_worker<update_state_mode>(&more);
      assert(more);

      switch (mode) {
      case update_state_mode::parameters: {
         state_id_ = get_manager()->messenger().receive_from_master_on_worker<RooFit::MultiProcess::State>(&more);
         assert(more);
         auto message = get_manager()->messenger().receive_from_master_on_worker<zmq::message_t>(&more);
         auto message_begin = message.data<update_state_t>();
         auto message_end = message_begin + message.size() / sizeof(update_state_t);
         std::vector<update_state_t> to_update(message_begin, message_end);
         for (auto const &item : to_update) {
            RooRealVar *rvar = (RooRealVar *)vars_.at(item.var_index);
            rvar->setVal(static_cast<double>(item.value));
            if (rvar->isConstant() != item.is_constant) {
               rvar->setConstant(static_cast<bool>(item.is_constant));
            }
         }

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

         break;
      }
      case update_state_mode::offsetting: {
         LikelihoodWrapper::enableOffsetting(get_manager()->messenger().receive_from_master_on_worker<bool>(&more));
         assert(!more);
         break;
      }
      }
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
            get_manager()->messenger().publish_from_master_to_workers(id_, update_state_mode::parameters, state_id_,
                                                                      std::move(message), std::move(offsets_message));
            offsets_previous_ = *component_offsets_;
         } else {
            get_manager()->messenger().publish_from_master_to_workers(id_, update_state_mode::parameters, state_id_,
                                                                      std::move(message));
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
      if (do_offset_ && component_offsets_->empty()) {
         likelihood_serial_->evaluate();
         // note: we don't need to get the offsets from the serial likelihood, because they are already coupled through
         // the shared_ptr
      }

      // update parameters that changed since last calculation (or creation if first time)
      updateWorkersParameters();

      // master fills queue with tasks
      auto N_tasks = getNEventTasks() * getNComponentTasks();
      for (std::size_t ix = 0; ix < N_tasks; ++ix) {
         get_manager()->queue()->add({id_, state_id_, ix});
      }
      n_tasks_at_workers_ = N_tasks;

      // wait for task results back from workers to master
      gather_worker_results();

      RooNaNPacker packedNaN;

      // Note: initializing result_ to results_[0] instead of zero-initializing it makes
      // a difference due to Kahan sum precision. This way, a single-worker run gives
      // the same result as a run with serial likelihood. Adding the terms to a zero
      // initial sum can cancel the carry in some cases, causing divergent values.
      result_ = results_[0];
      packedNaN.accumulate(results_[0].Sum());
      for (auto item_it = results_.cbegin() + 1; item_it != results_.cend(); ++item_it) {
         result_ += *item_it;
         packedNaN.accumulate(item_it->Sum());
      }
      results_.clear();

      if (packedNaN.getPayload() != 0) {
         result_ = ROOT::Math::KahanSum<double>(packedNaN.getNaNWithPayload());
      }

      if (TMath::IsNaN(result_.Sum())) {
         RooAbsReal::logEvalError(nullptr, GetName().c_str(), "function value is NAN");
      }
   }
}

// --- RESULT LOGISTICS ---

struct TaskEvalError {
   const RooAbsArg *originator = nullptr;
   std::string originatorName;
   std::string message;
   std::string serverValue; // a formatted string containing information about the originator's server object values

   TaskEvalError(const RooAbsArg *originator, std::string originatorName, std::string message, std::string serverValue)
      : originator(originator), originatorName(std::move(originatorName)), message(std::move(message)),
        serverValue(std::move(serverValue))
   {
   }

   explicit TaskEvalError(const zmq::message_t &zmqMessage)
   {
      auto data = reinterpret_cast<const unsigned char *>(zmqMessage.data());
      std::size_t dataIndex = 0;
      memcpy(&originator, &data[dataIndex], sizeof(originator));
      dataIndex += sizeof(originator);

      originatorName = decodeStringWithSize(data, dataIndex);
      message = decodeStringWithSize(data, dataIndex);
      serverValue = decodeStringWithSize(data, dataIndex);
   }

   explicit TaskEvalError(const unsigned char *zmqMessageData)
   {
      printf("in ctor zmqMessageData: %p\n", zmqMessageData);
      std::size_t dataIndex = 0;
      printf("hier\n");
      memcpy(&originator, &zmqMessageData[dataIndex], sizeof(originator));
      printf("ontvangen originator: %p\n", originator);
      dataIndex += sizeof(originator);
      printf("sizeof(originator): %zu\n", sizeof(originator));

      printf("ping\n");
      originatorName = decodeStringWithSize(zmqMessageData, dataIndex);
      printf("pong\n");
      message = decodeStringWithSize(zmqMessageData, dataIndex);
      printf("pang\n");
      serverValue = decodeStringWithSize(zmqMessageData, dataIndex);
      printf("pung\n");
   }

   static std::string decodeStringWithSize(const unsigned char *data, std::size_t &dataIndex)
   {
      std::size_t string_size = 0;
      printf("broink\n");
      memcpy(&string_size, &data[dataIndex], sizeof(string_size));
      dataIndex += sizeof(string_size);
      printf("string size: %zu\n", string_size);
      std::string decodedString{&reinterpret_cast<const char *>(data)[dataIndex], string_size};
      printf("blerv\n");
      dataIndex += string_size;
      return decodedString;
   }

   static void encodeStringWithSize(unsigned char *data, std::size_t &dataIndex, const std::string &string)
   {
      std::size_t string_size = string.size();
      printf("sizeof(string_size): %zu, sizeof(std::size_t): %zu\n", sizeof(string_size), sizeof(std::size_t));
      memcpy(&data[dataIndex], &string_size, sizeof(string_size));
      printf("dataIndex first: %zu\n", dataIndex);
      dataIndex += sizeof(string_size);
      printf("dataIndex then: %zu\n", dataIndex);
      memcpy(&data[dataIndex], string.data(), string.size());
      dataIndex += string.size();
      printf("dataIndex now: %zu\n", dataIndex);
   }

   std::size_t serializedSize() const
   {
      return sizeof(originator) + sizeof(originatorName.data()) * originatorName.size() +
             sizeof(message.data()) * message.size() + sizeof(serverValue.data()) * serverValue.size() +
             3 * sizeof(std::size_t); // the sizes of the strings must also be sent along
   }

   std::tuple<std::unique_ptr<unsigned char>, std::size_t> toZmqMessageBuffer() const
   {
      auto size = serializedSize();
      std::unique_ptr<unsigned char> data{reinterpret_cast<unsigned char *>(malloc(size))};
      std::size_t dataIndex = 0;

      printf("verstuurde originator: %p (sizeof(originator): %zu)\n", originator, sizeof(originator));
      memcpy(&data.get()[dataIndex], &originator, sizeof(originator));
      dataIndex += sizeof(originator);

      encodeStringWithSize(data.get(), dataIndex, originatorName);
      encodeStringWithSize(data.get(), dataIndex, message);
      encodeStringWithSize(data.get(), dataIndex, serverValue);
      printf("toZmqMessageBuffer pointer adres: %p\n", data.get());
      return {std::move(data), size};
   }

   static zmq::message_t vectorToZmqMessage(const std::vector<TaskEvalError> &evalErrors)
   {
      std::vector<std::tuple<std::unique_ptr<unsigned char>, std::size_t>> messageBuffers;
      std::size_t dataSize = sizeof(std::size_t); // we start out with the size of the vector
      for (auto &evalError : evalErrors) {
         messageBuffers.push_back(evalError.toZmqMessageBuffer());
         printf("toZmqMessageBuffer pointer adres in de messageBuffers vector: %p\n",
                std::get<0>(messageBuffers.back()).get());
         dataSize += std::get<1>(messageBuffers.back()) +
                     sizeof(std::size_t); // size of element block + one size indicator "header"
      }
      std::unique_ptr<unsigned char> data{reinterpret_cast<unsigned char *>(malloc(dataSize))};
      std::size_t vector_size = evalErrors.size();
      memcpy(&data.get()[0], &vector_size, sizeof(std::size_t));
      std::size_t dataIndex = sizeof(std::size_t);
      for (auto &messageBuffer : messageBuffers) {
         auto bufferSize = std::get<1>(messageBuffer);
         memcpy(&data.get()[dataIndex], &bufferSize, sizeof(std::size_t));
         dataIndex += sizeof(std::size_t);
         memcpy(&data.get()[dataIndex], std::get<0>(messageBuffer).get(), bufferSize);
         dataIndex += bufferSize;
      }
      assert(dataSize == dataIndex);
      printf("gemaakte zmqMessage op PID %d size: %zu\n", getpid(), dataSize);

      return {data.get(), dataSize};
   }

   static std::vector<TaskEvalError> zmqMessageToVector(const zmq::message_t &zmqMessage)
   {
      auto data = reinterpret_cast<const unsigned char *>(zmqMessage.data());
      printf("data ptr: %p\n", data);
      printf("data size: %zu\n", zmqMessage.size());
      std::size_t dataIndex = 0;
      std::size_t vector_size = 0;
      memcpy(&vector_size, &data[dataIndex], sizeof(vector_size));
      printf("vector size: %zu\n", vector_size);
      dataIndex += sizeof(vector_size);
      printf("data ptr + dataIndex: %p\n", data + dataIndex);

      std::vector<TaskEvalError> evalErrors;
      for (std::size_t element_ix = 0; element_ix < vector_size; ++element_ix) {
         printf("element %zu\n", element_ix);
         printf("in for loop: %p\n", &data[dataIndex]);

         evalErrors.emplace_back(&data[dataIndex]);
         printf("element %zu emplaced\n", element_ix);
         dataIndex += evalErrors.back().serializedSize();
         printf("element %zu klaar\n", element_ix);
         assert(dataIndex <= zmqMessage.size());
      }
      return evalErrors;
   }
};

void LikelihoodJob::send_back_task_result_from_worker(std::size_t /*task*/)
{
   int numErrors = RooAbsReal::numEvalErrors();
   std::vector<TaskEvalError> evalErrors;

   if (numErrors) {
//      printf("PID %d heeft errors, daar gaan we\n", getpid());
//      // Loop over errors
//      std::string objidstr;
//      {
//         std::ostringstream oss;
//         // Format string with object identity as this cannot be evaluated on the other side
//         oss << "LikelihoodJob on PID " << getpid() << "/";
//         objidstr = oss.str();
//      }
//      auto iter = RooAbsReal::evalErrorIter();
//      const RooAbsArg *ptr = nullptr;
//      for (int i = 0; i < RooAbsReal::numEvalErrorItems(); ++i) {
//         for (auto &item : iter->second.second) {
//            ptr = iter->first;
//            evalErrors.emplace_back(ptr, objidstr, item._msg, item._srvval);
//            printf("LikelihoodJob::send_back_task_result_from_worker(%s) sending error log Arg %p Msg %s\n",
//                   GetName().c_str(), iter->first, item._msg.c_str());
//         }
//      }
      // Clear error list on local side
      RooAbsReal::clearEvalErrorLog();
   }

   //   task_result_t task_result{id_, result_.Result(), result_.Carry(), !evalErrors.empty()};
   task_result_t task_result{id_, result_.Result(), result_.Carry(), numErrors > 0};
   zmq::message_t message(sizeof(task_result_t));
   memcpy(message.data(), &task_result, sizeof(task_result_t));

   //   if (evalErrors.empty()) {
//   printf("PID %d gaat beginnen met sturen zonder errors...\n", getpid());
   get_manager()->messenger().send_from_worker_to_master(std::move(message));
//   printf("...PID %d heeft gestuurd zonder errors\n", getpid());
   //
   //   } else {
   //      printf("PID %d gaat beginnen met sturen met errors...\n", getpid());
   //      zmq::message_t error_message = TaskEvalError::vectorToZmqMessage(evalErrors);
   //      printf("...PID %d heeft de error_message gemaakt...\n", getpid());
   //      get_manager()->messenger().send_from_worker_to_master(std::move(message), std::move(error_message));
   //      printf("...PID %d heeft gestuurd MET errors\n", getpid());
   //   }
}

bool LikelihoodJob::receive_task_result_on_master(const zmq::message_t &message)
{
//   printf("TESTSTSETSETSE\n");
//
//   std::string testString("pieter");
//
//   std::unique_ptr<unsigned char> testData{
//      reinterpret_cast<unsigned char *>(malloc(testString.size() + sizeof(std::size_t)))};
//
//   std::size_t testDataIndex = 0;
//   printf("testString size: %zu\n", testString.size());
//   TaskEvalError::encodeStringWithSize(testData.get(), testDataIndex, testString);
//   printf("testDataIndex: %zu\n", testDataIndex);
//   assert(testDataIndex == testString.size() + sizeof(std::size_t));
//   testDataIndex = 0;
//   auto testResult = TaskEvalError::decodeStringWithSize(testData.get(), testDataIndex);
//   printf("testDataIndex: %zu\n", testDataIndex);
//   assert(testDataIndex == testString.size() + sizeof(std::size_t));
//   printf("testResult: %s\n", testResult.c_str());
//   assert(testString == testResult);
//
//   printf("END TESTSTSETST\n");

   auto task_result = message.data<task_result_t>();
   results_.emplace_back(task_result->value, task_result->carry);
   if (task_result->has_errors) {
//      printf("we have some errors, let's goooo...\n");
      //      auto error_message = get_manager()->messenger().receive_from_worker_on_master<zmq::message_t>();
      //      printf("... got the message ...\n");
      //
      //      auto evalErrors = TaskEvalError::zmqMessageToVector(error_message);
      //      printf("... converted to vector ...\n");
      //      for (auto& evalError : evalErrors) {
      //         RooAbsReal::logEvalError(reinterpret_cast<const RooAbsReal *>(evalError.originator),
      //                                  evalError.originatorName.c_str(), evalError.message.c_str(),
      //                                  evalError.serverValue.c_str());
      //      }
      //      printf("... put them all back into the error log!\n");
//      printf("...let's just put one error to rule them all in the log\n");
      RooAbsReal::logEvalError(nullptr, "LikelihoodJob", "evaluation errors at the worker processes", "no servervalue");
   }
   --n_tasks_at_workers_;
   bool job_completed = (n_tasks_at_workers_ == 0);
   return job_completed;
}

// --- END OF RESULT LOGISTICS ---

void LikelihoodJob::evaluate_task(std::size_t task)
{
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
      RooNaNPacker packedNaN;
      for (std::size_t comp_ix = components_first; comp_ix < components_last; ++comp_ix) {
         auto component_result = likelihood_->evaluatePartition({section_first, section_last}, comp_ix, comp_ix + 1);
         packedNaN.accumulate(component_result.Sum());
         if (do_offset_ && section_last == 1 && (*component_offsets_)[comp_ix] != ROOT::Math::KahanSum<double>(0, 0)) {
            // we only subtract at the end of event sections, otherwise the offset is subtracted for each event split
            result_ += (component_result - (*component_offsets_)[comp_ix]);
         } else {
            result_ += component_result;
         }
      }
      if (packedNaN.getPayload() != 0) {
         result_ = ROOT::Math::KahanSum<double>(packedNaN.getNaNWithPayload());
      }

      break;
   }
   }
}

void LikelihoodJob::enableOffsetting(bool flag)
{
   likelihood_serial_->enableOffsetting(flag);
   LikelihoodWrapper::enableOffsetting(flag);
   if (RooFit::MultiProcess::JobManager::is_instantiated()) {
      printf("WARNING: when calling MinuitFcnGrad::setOffsetting after the run has already been started the "
             "MinuitFcnGrad::likelihood_in_gradient object (a LikelihoodSerial) on the workers can no longer be "
             "updated! This function (LikelihoodJob::enableOffsetting) can in principle be used outside of "
             "MinuitFcnGrad, but be aware of this limitation. To do a minimization with a different offsetting "
             "setting, please delete all RooFit::MultiProcess based objects so that the forked processes are killed "
             "and then set up a new RooMinimizer.\n");
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
