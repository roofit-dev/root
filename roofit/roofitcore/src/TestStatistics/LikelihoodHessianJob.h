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

#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodHessianJob
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodHessianJob

#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/MinimumState.h"

#include "RooFit/MultiProcess/Job.h"

namespace RooFit {
namespace TestStatistics {
ROOT::Minuit2::MinimumState hessian_calculator(const ROOT::Minuit2::MnStrategy &strategy,
		                               const ROOT::Minuit2::MnFcn &mfcn,
					       const ROOT::Minuit2::MinimumState &st,
                                               const ROOT::Minuit2::MnUserTransformation &trafo,
					       unsigned int maxcalls);
} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_TESTSTATISTICS_LikelihoodHessianJob

