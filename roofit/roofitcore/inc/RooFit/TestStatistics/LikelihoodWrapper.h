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

#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodWrapper
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodWrapper

#include "RooArgSet.h"
#include "RooAbsArg.h" // enum ConstOpCode

#include <Fit/ParameterSettings.h>
#include <Math/MinimizerOptions.h>
#include <Math/Util.h>

#include <memory> // shared_ptr
#include <string>

// forward declaration
class RooMinimizer;

namespace RooFit {
namespace TestStatistics {

// forward declaration
class RooAbsL;

/// For communication with wrappers, an instance of this struct must be shared between them and MinuitFcnGrad. It keeps
/// track of what has been evaluated for the current parameter set provided by Minuit.
struct WrapperCalculationCleanFlags {
   // indicate whether that part has been calculated since the last parameter update
   bool likelihood = false;
   bool gradient = false;

   void set_all(bool value)
   {
      likelihood = value;
      gradient = value;
   }
};

enum class LikelihoodType { unbinned, binned, subsidiary, sum };

enum class LikelihoodMode { serial, multiprocess };

/// Previously, offsetting was only implemented for RooNLLVar components of a likelihood,
/// not for RooConstraintSum terms. To emulate this behavior, use OffsettingMode::legacy. To
/// also offset the RooSubsidiaryL component (equivalent of RooConstraintSum) of RooSumL
/// likelihoods, use OffsettingMode::full.
enum class OffsettingMode { legacy, full };

class LikelihoodWrapper {
public:
   LikelihoodWrapper(std::shared_ptr<RooAbsL> likelihood,
                     std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean);
   LikelihoodWrapper(std::shared_ptr<RooAbsL> likelihood,
                     std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean,
                     std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> offsets,
                     std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> offsets_save);
   virtual ~LikelihoodWrapper() = default;
   virtual LikelihoodWrapper *clone() const = 0;

   static std::unique_ptr<LikelihoodWrapper> create(LikelihoodMode likelihoodMode, std::shared_ptr<RooAbsL> likelihood,
                                                    std::shared_ptr<WrapperCalculationCleanFlags> calculationIsClean);
   static std::unique_ptr<LikelihoodWrapper> create(LikelihoodMode likelihoodMode, std::shared_ptr<RooAbsL> likelihood,
                                                    std::shared_ptr<WrapperCalculationCleanFlags> calculationIsClean,
                                                    std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> offsets,
                                                    std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> offsets_save);

   /// \brief Triggers (possibly asynchronous) evaluation of the likelihood
   ///
   /// In parallel strategies, it may be advantageous to allow a calling process to continue on with other tasks while
   /// the calculation is offloaded to another process or device, like a GPU. For this reason, evaluate() does not
   /// return the result, this is done in getResult().
   virtual void evaluate() = 0;
   /// \brief Return the latest result of a likelihood evaluation.
   ///
   /// Returns the result that was stored after calling evaluate(). It is up to the implementer to make sure the stored
   /// value represents the most recent evaluation call, e.g. by using a mutex.
   virtual ROOT::Math::KahanSum<double> getResult() const = 0;

   /// Synchronize minimizer settings with calculators in child classes
   virtual void synchronizeWithMinimizer(const ROOT::Math::MinimizerOptions &options);
   virtual void synchronizeParameterSettings(const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings);
   /// Minuit passes in parameter values that may not conform to RooFit internal standards (like applying range
   /// clipping), but that the specific calculator does need. This function can be implemented to receive these
   /// Minuit-internal values:
   virtual void updateMinuitInternalParameterValues(const std::vector<double> &minuit_internal_x);
   virtual void updateMinuitExternalParameterValues(const std::vector<double> &minuit_external_x);

   // The following functions are necessary from MinuitFcnGrad to reach likelihood properties:
   void constOptimizeTestStatistic(RooAbsArg::ConstOpCode opcode, bool doAlsoTrackingOpt);
   double defaultErrorLevel() const;
   virtual std::string GetName() const;
   virtual std::string GetTitle() const;
   inline virtual bool isOffsetting() const { return do_offset_; }
   virtual void enableOffsetting(bool flag);
   void setOffsettingMode(OffsettingMode mode);
   inline std::vector<ROOT::Math::KahanSum<double>> offsets() const { return *component_offsets_; }
   void setApplyWeightSquared(bool flag);

   void setLogfile(std::shared_ptr<std::ofstream> logfile);
   void initializeLogfile(std::shared_ptr<std::ofstream>& logfile);
protected:
   std::shared_ptr<RooAbsL> likelihood_;
   LikelihoodType likelihood_type_;
   std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean_;

   bool do_offset_ = false;
   std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> component_offsets_;
   std::shared_ptr<std::vector<ROOT::Math::KahanSum<double>>> component_offsets_save_;
   void calculate_offsets();
   void clearOffsets();
   OffsettingMode offsetting_mode_ = OffsettingMode::legacy;
   void swapOffsets(const std::vector<std::size_t>& component_indices);

   std::shared_ptr<std::ofstream> logfile_;
};

} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_TESTSTATISTICS_LikelihoodWrapper
