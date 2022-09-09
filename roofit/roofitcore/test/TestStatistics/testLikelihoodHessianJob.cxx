/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2022, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include "RooRandom.h"
#include "RooWorkspace.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooDataHist.h" // complete type in Binned test
#include "RooCategory.h" // complete type in MultiBinnedConstraint test
#include "RooFit/TestStatistics/RooUnbinnedL.h"
#include "RooFit/TestStatistics/RooBinnedL.h"
#include "RooFit/TestStatistics/optional_parameter_types.h"
#include "RooFit/TestStatistics/buildLikelihood.h"
#include "RooFit/TestStatistics/RooRealL.h"
#include "RooFit/MultiProcess/Config.h"
#include "../../src/RooGradMinimizerFcn.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FCNAdapter.h"
#include "Minuit2/FCNGradAdapter.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnFcn.h"  // needed to complete type for lambda to MnHesse::calculator
#include "Math/Util.h" // KahanSum
#include "TMatrixDSym.h"

#include "../../src/TestStatistics/LikelihoodHessianJob.h"

#include <stdexcept> // runtime_error

#include "gtest/gtest.h"
#include "../test_lib.h" // generate_1D_gaussian_pdf_nll

class Environment : public testing::Environment {
public:
   void SetUp() override
   {
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
      RooFit::MultiProcess::Config::setDefaultNWorkers(2);
      ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   }
};

// Previously, we just called AddGlobalTestEnvironment in global namespace, but this caused either a warning about an
// unused declared variable (the return value of the call) or a parsing problem that the compiler can't handle if you
// don't store the return value at all. The solution is to just define this manual main function. The default gtest
// main function does InitGoogleTest and RUN_ALL_TESTS, we add the environment call in the middle.
int main(int argc, char **argv)
{
   testing::InitGoogleTest(&argc, argv);
   testing::AddGlobalTestEnvironment(new Environment);
   return RUN_ALL_TESTS();
}

/// Create an IMultiGradFunction from a FCNGradAdapter, assuming it is a likelihood
ROOT::Minuit2::FCNGradAdapter<ROOT::Math::IMultiGradFunction> to_adapter(const ROOT::Math::IMultiGradFunction &func)
{
   // errorDef is typically 0.5 for likelihoods
   auto errorDef = 0.5;
   return ROOT::Minuit2::FCNGradAdapter<ROOT::Math::IMultiGradFunction>(func, errorDef);
}

void printMatrix(const ROOT::Minuit2::LASymMatrix &matrix)
{
   for (std::size_t i_row = 0; i_row < matrix.Nrow(); ++i_row) {
      for (std::size_t i_col = 0; i_col < matrix.Ncol(); ++i_col) {
         printf("%f ", matrix(i_row, i_col));
      }
      printf("\n");
   }
}

void printMatrix(const TMatrixDSym &matrix)
{
   for (std::size_t i_row = 0; i_row < static_cast<std::size_t>(matrix.GetNrows()); ++i_row) {
      for (std::size_t i_col = 0; i_col < static_cast<std::size_t>(matrix.GetNcols()); ++i_col) {
         printf("%f ", matrix(i_row, i_col));
      }
      printf("\n");
   }
}


// Due to the integrated nature of Minuit2, it is very hard to use MnHesse
// outside of an actual Minuit2 minimization run. This suite therefore simply
// runs two minimizations using RooMinimizer, switching out the hessian
// calculator function between the two runs and checking the results to make
// sure our custom hessian calculator works.
class MnHesseCalculatorTest: public ::testing::Test {
protected:
   void SetUp() override
   {
      RooRandom::randomGenerator()->SetSeed(seed);

      // parameters
      unsigned int N = 2;

      std::tie(nll, pdf, data, values) = generate_ND_gaussian_pdf_nll(w, N, 1000);

      RooArgSet *savedValues = dynamic_cast<RooArgSet *>(values->snapshot());
      if (savedValues == nullptr) {
         throw std::runtime_error("params->snapshot() cannot be casted to RooArgSet!");
      }

      // --------

      minimizer = std::make_unique<RooMinimizer>(*nll, RooMinimizer::FcnMode::gradient);

      minimizer->setStrategy(0);
      minimizer->setPrintLevel(-1);

      // first minimize to get to a point where the Hessian will be sane
      minimizer->migrad();

      // then do the Hessian
      minimizer->hesse();

      RooFitResult *m0result = minimizer->lastMinuitFit();
      cov0 = std::make_unique<TMatrixDSym>(m0result->covarianceMatrix());

      // --------

      *values = *savedValues;
   }

   void TearDown() override
   {
      // switch back to default_calculator
      ROOT::Minuit2::MnHesse::calculator = ROOT::Minuit2::MnHesse::default_calculator;
   }

   std::size_t seed = 39;
   RooWorkspace w;
   std::unique_ptr<RooAbsReal> nll;
   std::unique_ptr<RooArgSet> values;
   RooAbsPdf *pdf;
   RooAbsData *data;
   std::shared_ptr<RooFit::TestStatistics::RooAbsL> likelihood;
   std::unique_ptr<RooMinimizer> minimizer;
   std::unique_ptr<TMatrixDSym> cov0;
};

TEST_F(MnHesseCalculatorTest, DummyLambda)
{
   // In this test we just switch the calculator out with a dummy lambda that
   // prints a line and then forwards to the default calculator.

   ROOT::Minuit2::MnHesse::calculator = [](auto&& ...args) {
      printf("hi there, going to run the default_calculator now from this lambda\n");
      return ROOT::Minuit2::MnHesse::default_calculator(std::forward<decltype(args)>(args)...);
   };

   minimizer->hesse();

   RooFitResult *m1result = minimizer->lastMinuitFit();
   auto cov1 = m1result->covarianceMatrix();

   EXPECT_EQ(*cov0, cov1);
}

TEST_F(MnHesseCalculatorTest, MPCalculator)
{
   // Here we switch the calculator out with the MultiProcess parallelized
   // version of the default calculator.

   ROOT::Minuit2::MnHesse::calculator = [](auto&& ...args) {
      printf("hi there, going to run the RooFit::TestStatistics::hessian_calculator now from this lambda\n");
      return RooFit::TestStatistics::hessian_calculator(std::forward<decltype(args)>(args)...);
   };

   minimizer->hesse();

   RooFitResult *m1result = minimizer->lastMinuitFit();
   auto cov1 = m1result->covarianceMatrix();

   EXPECT_EQ(*cov0, cov1);
}

