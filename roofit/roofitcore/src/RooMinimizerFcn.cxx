/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl   *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
/// \class RooMinimizerFcn
/// RooMinimizerFcn is an interface to the ROOT::Math::IBaseFunctionMultiDim,
/// a function that ROOT's minimisers use to carry out minimisations.
///

#include "RooMinimizerFcn.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooNaNPacker.h"

#include "Math/Functor.h"
#include "TMatrixDSym.h"
#include "RooAddition.h"

#include <fstream>
#include <iomanip>

using namespace std;

namespace {

// Helper function that wraps RooAbsArg::getParameters and directly returns the
// output RooArgSet. To be used in the initializer list of the RooMinimizerFcn
// constructor.
RooArgSet getParameters(RooAbsReal const &funct)
{
   RooArgSet out;
   funct.getParameters(nullptr, out);
   return out;
}

} // namespace

// use reference wrapper for the Functor, such that the functor points to this RooMinimizerFcn by reference.
RooMinimizerFcn::RooMinimizerFcn(RooAbsReal *funct, RooMinimizer *context)
   : RooAbsMinimizerFcn(getParameters(*funct), context), _funct(funct)
{
    if (context->_cfg.useGradient && funct->hasGradient()) {
        _multiGenFcn = std::make_unique<ROOT::Math::GradFunctor>(this, &RooMinimizerFcn::operator(),
                                                                 &RooMinimizerFcn::evaluateGradient, getNDim());
    } else {
        _multiGenFcn = std::make_unique<ROOT::Math::Functor>(std::cref(*this), getNDim());
    }

   do_eval_log_file_ = std::make_shared<std::ofstream>("RooMinimizerFcnDoEvalLogFile.csv");
   *do_eval_log_file_ << std::setprecision(17);
   for (auto &parameter : context->fitter()->Config().ParamsSettings()) {
      *do_eval_log_file_ << parameter.Name() << "," << std::flush;
   }

   // we also need component names in the log file for simultaneous pdf fits:
   auto possibleNLL = dynamic_cast<RooAddition*>((*funct->ownedComponents())[0]);
   if (possibleNLL != nullptr) {
      for (auto& component : possibleNLL->list()) {
         *do_eval_log_file_ << component->GetName() << "," << std::flush;
      }
   }

   *do_eval_log_file_ << "fcn_value" << std::endl;
}

ROOT::Math::IBaseFunctionMultiDim *RooMinimizerFcn::Clone() const
{
   return new RooMinimizerFcn(*this);
}

void RooMinimizerFcn::setOptimizeConstOnFunction(RooAbsArg::ConstOpCode opcode, bool doAlsoTrackingOpt)
{
   _funct->constOptimizeTestStatistic(opcode, doAlsoTrackingOpt);
}

/// Evaluate function given the parameters in `x`.
double RooMinimizerFcn::operator()(const double *x) const
{
   for (std::size_t ix = 0; ix < NDim(); ++ix) {
      *do_eval_log_file_ << x[ix] << "," << std::flush;
   }

   // Set the parameter values for this iteration
   for (unsigned index = 0; index < _nDim; index++) {
      if (_logfile)
         (*_logfile) << x[index] << " ";
      SetPdfParamVal(index, x[index]);
   }

   // Calculate the function for these parameters
   RooAbsReal::setHideOffset(false);
   double fvalue = _funct->getVal();

   ///// BEGIN temporary for logging:
   auto possibleNLL = dynamic_cast<RooAddition*>((*_funct->ownedComponents())[0]);
   if (possibleNLL != nullptr) {
      auto& set = possibleNLL->list();
      const RooArgSet* nset = set.nset();

      for (auto& component : set) {
         *do_eval_log_file_ << component->getVal(nset) << "," << std::flush;
      }
   }
   ///// END temporary for logging

   RooAbsReal::setHideOffset(true);

   if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {
      printEvalErrors();
      RooAbsReal::clearEvalErrorLog();
      _numBadNLL++;

      if (cfg().doEEWall) {
         const double badness = RooNaNPacker::unpackNaN(fvalue);
         fvalue = (std::isfinite(_maxFCN) ? _maxFCN : 0.) + cfg().recoverFromNaN * badness;
      }
   } else {
      if (_evalCounter > 0 && _evalCounter == _numBadNLL) {
         // This is the first time we get a valid function value; while before, the
         // function was always invalid. For invalid  cases, we returned values > 0.
         // Now, we offset valid values such that they are < 0.
         _funcOffset = -fvalue;
      }
      fvalue += _funcOffset;
      _maxFCN = std::max(fvalue, _maxFCN);
   }

   // Optional logging
   if (_logfile)
      (*_logfile) << setprecision(15) << fvalue << setprecision(4) << endl;
   if (cfg().verbose) {
      cout << "\nprevFCN" << (_funct->isOffsetting() ? "-offset" : "") << " = " << setprecision(10) << fvalue
           << setprecision(4) << "  ";
      cout.flush();
   }

   finishDoEval();

   *do_eval_log_file_ << fvalue << std::endl;
   return fvalue;
}

void RooMinimizerFcn::evaluateGradient(const double *x, double *out) const
{
   // Set the parameter values for this iteration
   for (unsigned index = 0; index < _nDim; index++) {
      if (_logfile)
         (*_logfile) << x[index] << " ";
      SetPdfParamVal(index, x[index]);
   }

   _funct->gradient(out);

   // Optional logging
   if (cfg().verbose) {
      std::cout << "\n    gradient = ";
      for (std::size_t i = 0; i < getNDim(); ++i) {
         std::cout << out[i] << ", ";
      }
   }
}

std::string RooMinimizerFcn::getFunctionName() const
{
   return _funct->GetName();
}

std::string RooMinimizerFcn::getFunctionTitle() const
{
   return _funct->GetTitle();
}

void RooMinimizerFcn::setOffsetting(bool flag)
{
   _funct->enableOffsetting(flag);
}
