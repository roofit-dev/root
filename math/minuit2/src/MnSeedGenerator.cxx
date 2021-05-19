// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 * Copyright (c) 2017 Patrick Bos, Netherlands eScience Center        *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnSeedGenerator.h"
#include "Minuit2/MinimumSeed.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/GradientCalculator.h"
#include "Minuit2/InitialGradientCalculator.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/FunctionGradient.h"
#include "Minuit2/MinimumError.h"
#include "Minuit2/MnMatrix.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MinuitParameter.h"
#include "Minuit2/MnLineSearch.h"
#include "Minuit2/MnParabolaPoint.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/VariableMetricEDMEstimator.h"
#include "Minuit2/NegativeG2LineSearch.h"
#include "Minuit2/AnalyticalGradientCalculator.h"
#include "Minuit2/Numerical2PGradientCalculator.h"
#include "Minuit2/HessianGradientCalculator.h"
#include "Minuit2/MnPrint.h"

#include <cmath>


namespace ROOT {

   namespace Minuit2 {


MinimumSeed MnSeedGenerator::operator()(const MnFcn& fcn, const GradientCalculator& gc, const MnUserParameterState& st, const MnStrategy& stra) const {

   MnPrint print("MnSeedGenerator");

   // find seed (initial minimization point) using the calculated gradient
   const unsigned int n = st.VariableParameters();
   const MnMachinePrecision& prec = st.Precision();

   print.Debug(n, "free parameters, FCN pointer", &fcn);

   // initial starting values
   MnAlgebraicVector x(n);
   for(unsigned int i = 0; i < n; i++) x(i) = st.IntParameters()[i];

//  std::cout << "-- MnSeedGenerator::operator(.., GradientCalculator, ..):" << std::endl;
//  for (int i = 0; i < n; ++i) {
//    std::cout << std::hexfloat << "x=("<< x(i) << ",\t";
//  }
//  std::cout << ")" << std::endl;

  double fcnmin = fcn(x);

   MinimumParameters pa(x, fcnmin);
//  std::cout << "... doing gc(pa) ..." << std::endl;
//  std::cout << "-- hier? 1 --" << std::endl;
   FunctionGradient dgrad = gc(pa);
//  std::cout << "-- hier? 2 --" << std::endl;
//  std::cout << "dgrad.Vec: " << dgrad.Vec() << std::endl;
//  std::cout << "dgrad.G2: " << dgrad.G2() << std::endl;
   MnAlgebraicSymMatrix mat(n);
   double dcovar = 1.;
   if(st.HasCovariance()) {
//     std::cout << "has covariance" << std::endl;
      for(unsigned int i = 0; i < n; i++)
         for(unsigned int j = i; j < n; j++) mat(i,j) = st.IntCovariance()(i,j);
      dcovar = 0.;
   } else {
//     std::cout << "has no covariance" << std::endl;
      for(unsigned int i = 0; i < n; i++)
         mat(i,i) = (std::fabs(dgrad.G2()(i)) > prec.Eps2() ? 1./dgrad.G2()(i) : 1.);
   }
   MinimumError err(mat, dcovar);

   double edm = VariableMetricEDMEstimator().Estimate(dgrad, err);
   MinimumState state(pa, err, dgrad, edm, fcn.NumOfCalls());

   print.Info("Initial state:", MnPrint::Oneline(state));

   NegativeG2LineSearch ng2ls;
   if(ng2ls.HasNegativeG2(dgrad, prec)) {
      print.Debug("Negative G2 Found",
        "\n  point:", x,
        "\n  grad :", dgrad.Grad(),
        "\n  g2   :", dgrad.G2());

      state = ng2ls(fcn, state, gc, prec);

      print.Info("Negative G2 found - new state:", state);
   }


   if(stra.Strategy() == 2 && !st.HasCovariance()) {
      //calculate full 2nd derivative

      print.Debug("calling MnHesse");

      MinimumState tmp = MnHesse(stra)(fcn, state, st.Trafo());

      print.Info("run Hesse - new state:", tmp);

//     std::cout << "-- MnSeedGenerator::operator(.., GradientCalculator, ..), strategy 2 end state:" << std::endl;
//     for (int i = 0; i < n; ++i) {
//       std::cout << std::hexfloat << "x=("<< tmp.Vec()(i) << ",\t";
//     }
//     std::cout << ")" << std::endl;

     return MinimumSeed(tmp, st.Trafo());
   }

//  std::cout << "-- MnSeedGenerator::operator(.., GradientCalculator, ..), regular end state:" << std::endl;
//  for (int i = 0; i < n; ++i) {
//    std::cout << std::hexfloat << "x=("<< state.Vec()(i) << ",\t";
//  }
//  std::cout << ")" << std::endl;

  return MinimumSeed(state, st.Trafo());
}


MinimumSeed MnSeedGenerator::operator()(const MnFcn& fcn, const AnalyticalGradientCalculator& gc, const MnUserParameterState& st, const MnStrategy& stra) const {
   MnPrint print("MnSeedGenerator");

   // find seed (initial point for minimization) using analytical gradient
   unsigned int n = st.VariableParameters();
   const MnMachinePrecision& prec = st.Precision();

   int printLevel = MnPrint::Level();

  // initial starting values
   MnAlgebraicVector x(n);
   for(unsigned int i = 0; i < n; i++) x(i) = st.IntParameters()[i];

//  std::cout << "-- MnSeedGenerator::operator(.., AnalyticalGradientCalculator, ..):" << std::endl;
//  for (int i = 0; i < n; ++i) {
//    std::cout << std::hexfloat << "x=("<< x(i) << ",\t";
//  }
//  std::cout << ")" << std::endl;

  double fcnmin = fcn(x);
   MinimumParameters pa(x, fcnmin);

//  std::cout << "... creating igc ..." << std::endl;
//   InitialGradientCalculator igc(fcn, st.Trafo(), stra);
//  std::cout << "... doing igc(pa) ..." << std::endl;
//   FunctionGradient tmp = igc(pa);
//  std::cout << "tmp.G2: " << tmp.G2() << std::endl;
//  std::cout << "-- hier? 1 --" << std::endl;
   FunctionGradient grd = gc(pa);
//  std::cout << "grd.G2: " << grd.G2() << std::endl;
//  std::cout << "-- hier? 2 --" << std::endl;

//  FunctionGradient dgrad(grd.Grad(), tmp.G2(), tmp.Gstep());
    FunctionGradient dgrad(grd.Grad(), grd.G2(), grd.Gstep());
//  std::cout << "dgrad.Vec: " << dgrad.Vec() << std::endl;
//  std::cout << "dgrad.G2: " << dgrad.G2() << std::endl;

   if(gc.CheckGradient()) {
      bool good = true;
      HessianGradientCalculator hgc(fcn, st.Trafo(), MnStrategy(2));
      std::pair<FunctionGradient, MnAlgebraicVector> hgrd = hgc.DeltaGradient(pa, dgrad);
      for(unsigned int i = 0; i < n; i++) {
         if(std::fabs(hgrd.first.Grad()(i) - grd.Grad()(i)) > hgrd.second(i)) {
            int externalParameterIndex = st.Trafo().ExtOfInt(i);
            const char * parameter_name = st.Trafo().Name(externalParameterIndex);
            print.Warn("Gradient discrepancy of external Parameter too large:"
              "parameter_name =", parameter_name,
              "externalParameterIndex =", externalParameterIndex,
              "internal =", i
            );
            good = false;
         }
      }
      if(!good) {
         print.Error("Minuit does not accept user specified Gradient. To force acceptance, override 'virtual bool CheckGradient() const' of FCNGradientBase.h in the derived class.");

         assert(good);
      }
   }

   MnAlgebraicSymMatrix mat(n);
   double dcovar = 1.;
   if(st.HasCovariance()) {
//     std::cout << "has covariance" << std::endl;
      for(unsigned int i = 0; i < n; i++)
         for(unsigned int j = i; j < n; j++) mat(i,j) = st.IntCovariance()(i,j);
      dcovar = 0.;
   } else {
//     std::cout << "has no covariance" << std::endl;
      for(unsigned int i = 0; i < n; i++)
         mat(i,i) = (std::fabs(dgrad.G2()(i)) > prec.Eps2() ? 1./dgrad.G2()(i) : 1.);
   }
   MinimumError err(mat, dcovar);
   double edm = VariableMetricEDMEstimator().Estimate(dgrad, err);
   MinimumState state(pa, err, dgrad, edm, fcn.NumOfCalls());

   NegativeG2LineSearch ng2ls;
//  std::cout << dgrad.Vec() << std::endl;
//  std::cout << prec << std::endl;
//  std::cout << ng2ls.HasNegativeG2(dgrad, prec) << std::endl;
//  std::cout << dgrad.Vec().size() << std::endl;

   if(ng2ls.HasNegativeG2(dgrad, prec)) {
//      Numerical2PGradientCalculator ngc(fcn, st.Trafo(), stra);
//      state = ng2ls(fcn, state, ngc, prec);
      state = ng2ls(fcn, state, gc, prec);

     if (printLevel >1) {
       MnPrint::PrintState(std::cout, state, "MnSeedGenerator: Negative G2 found - new state:  ");
     }

   }

   if(stra.Strategy() == 2 && !st.HasCovariance()) {
      //calculate full 2nd derivative
      MinimumState tmpState = MnHesse(stra)(fcn, state, st.Trafo());

//     std::cout << "-- MnSeedGenerator::operator(.., AnalyticalGradientCalculator, ..), strategy 2 end state:" << std::endl;
//     for (int i = 0; i < n; ++i) {
//       std::cout << std::hexfloat << "x=("<< tmpState.Vec()(i) << ",\t";
//     }
//     std::cout << ")" << std::endl;

     return MinimumSeed(tmpState, st.Trafo());
   }

//  std::cout << "-- MnSeedGenerator::operator(.., AnalyticalGradientCalculator, ..), regular end state:" << std::endl;
//  for (int i = 0; i < n; ++i) {
//    std::cout << std::hexfloat << "x=("<< state.Vec()(i) << ",\t";
//  }
//  std::cout << ")" << std::endl;

  return MinimumSeed(state, st.Trafo());
}

   }  // namespace Minuit2

}  // namespace ROOT
