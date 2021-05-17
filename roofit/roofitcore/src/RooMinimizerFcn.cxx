/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOMINIMIZER

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
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooGaussMinimizer.h"

#include "TClass.h"
#include "TMatrixDSym.h"

#include <fstream>
#include <iomanip>

using namespace std;

RooMinimizerFcn::RooMinimizerFcn(RooAbsReal *funct, RooMinimizer* context,
			   bool verbose) :
  RooAbsMinimizerFcn(*funct->getParameters(RooArgSet()), context, verbose), _funct(funct)
{}



RooMinimizerFcn::RooMinimizerFcn(const RooMinimizerFcn& other) : RooAbsMinimizerFcn(other), ROOT::Math::IBaseFunctionMultiDim(other),
  _funct(other._funct)
{}


RooMinimizerFcn::~RooMinimizerFcn()
{}


ROOT::Math::IBaseFunctionMultiDim* RooMinimizerFcn::Clone() const 
{  
  return new RooMinimizerFcn(*this) ;
}


void RooMinimizerFcn::optimizeConstantTerms(bool constStatChange, bool constValChange) {
   if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;

      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: set of constant parameters changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ConfigChange) ;
   } else if (constValChange) {
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: constant parameter values changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ValueChange) ;
   }

   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
}


/// Evaluate function given the parameters in `x`.
double RooMinimizerFcn::DoEval(const double *x) const {

  // Set the parameter values for this iteration
  for (int index = 0; index < _nDim; index++) {
    if (_logfile) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // Calculate the function for these parameters  
  RooAbsReal::setHideOffset(kFALSE) ;
  double fvalue = _funct->getVal();
  RooAbsReal::setHideOffset(kTRUE) ;

  if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {
    printEvalErrors();
    RooAbsReal::clearEvalErrorLog() ;
    _numBadNLL++ ;

    if (_doEvalErrorWall) {
      const double badness = RooNaNPacker::unpackNaN(fvalue);
      fvalue = (std::isfinite(_maxFCN) ? _maxFCN : 0.) + _recoverFromNaNStrength * badness;
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
  if (_verbose) {
    cout << "\nprevFCN" << (_funct->isOffsetting()?"-offset":"") << " = " << setprecision(10) 
         << fvalue << setprecision(4) << "  " ;
    cout.flush() ;
  }

  _evalCounter++ ;

  return fvalue;
}

#endif

