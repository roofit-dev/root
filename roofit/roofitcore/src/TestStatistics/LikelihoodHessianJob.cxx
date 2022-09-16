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

#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/Numerical2PGradientCalculator.h"
#include "Minuit2/HessianGradientCalculator.h"
#include "Minuit2/MnPosDef.h"
#include "Minuit2/VariableMetricEDMEstimator.h"

#include "LikelihoodHessianJob.h"


namespace RooFit {
namespace TestStatistics {

using ROOT::Minuit2::MinimumState;
using ROOT::Minuit2::MnStrategy;
using ROOT::Minuit2::MnFcn;
using ROOT::Minuit2::MnUserTransformation;
using ROOT::Minuit2::MnMachinePrecision;
using ROOT::Minuit2::MinimumError;
using ROOT::Minuit2::FunctionGradient;
using ROOT::Minuit2::MnAlgebraicVector;
using ROOT::Minuit2::MnAlgebraicSymMatrix;
using ROOT::Minuit2::MnPrint;
using ROOT::Minuit2::HessianGradientCalculator;
using ROOT::Minuit2::Numerical2PGradientCalculator;
using ROOT::Minuit2::MnPosDef;
using ROOT::Minuit2::VariableMetricEDMEstimator;

/// Calculates the hessian from Minuit2 using LikelihoodHessianJob
///
/// This function can be used to replace MnHesse::calculator. Minuit2 will
/// then use it internally during minimization and during explicit calls to
/// calculate the hessian, e.g. from RooMinimizer::hesse().
MinimumState hessian_calculator(const MnStrategy &strategy, const MnFcn &mfcn, const MinimumState &st,
                                const MnUserTransformation &trafo, unsigned int maxcalls)
{
   MnPrint print("RooFit::MultiProcess enabled MnHesse calculator replacement");

   const MnMachinePrecision &prec = trafo.Precision();
   // make sure starting at the right place
   double amin = mfcn(st.Vec());
   double aimsag = std::sqrt(prec.Eps2()) * (std::fabs(amin) + mfcn.Up());

   // diagonal Elements first

   unsigned int n = st.Parameters().Vec().size();
   if (maxcalls == 0)
      maxcalls = 200 + 100 * n + 5 * n * n;

   MnAlgebraicSymMatrix vhmat(n);
   MnAlgebraicVector g2 = st.Gradient().G2();
   MnAlgebraicVector gst = st.Gradient().Gstep();
   MnAlgebraicVector grd = st.Gradient().Grad();
   MnAlgebraicVector dirin = st.Gradient().Gstep();
   MnAlgebraicVector yy(n);

   // case gradient is not numeric (could be analytical or from FumiliGradientCalculator)

   if (st.Gradient().IsAnalytical()) {
      Numerical2PGradientCalculator igc(mfcn, trafo, strategy);
      FunctionGradient tmp = igc(st.Parameters());
      gst = tmp.Gstep();
      dirin = tmp.Gstep();
      g2 = tmp.G2();
   }

   MnAlgebraicVector x = st.Parameters().Vec();

   print.Debug("Gradient is", st.Gradient().IsAnalytical() ? "analytical" : "numerical", "\n  point:", x,
               "\n  fcn  :", amin, "\n  grad :", grd, "\n  step :", gst, "\n  g2   :", g2);

   for (unsigned int i = 0; i < n; i++) {

      double xtf = x(i);
      double dmin = 8. * prec.Eps2() * (std::fabs(xtf) + prec.Eps2());
      double d = std::fabs(gst(i));
      if (d < dmin)
         d = dmin;

      print.Debug("Derivative parameter", i, "d =", d, "dmin =", dmin);

      for (unsigned int icyc = 0; icyc < strategy.HessianNCycles(); icyc++) {
         double sag = 0.;
         double fs1 = 0.;
         double fs2 = 0.;
         for (unsigned int multpy = 0; multpy < 5; multpy++) {
            x(i) = xtf + d;
            fs1 = mfcn(x);
            x(i) = xtf - d;
            fs2 = mfcn(x);
            x(i) = xtf;
            sag = 0.5 * (fs1 + fs2 - 2. * amin);

            print.Debug("cycle", icyc, "mul", multpy, "\tsag =", sag, "d =", d);

            //  Now as F77 Minuit - check that sag is not zero
            if (sag != 0)
               goto L30; // break
            if (trafo.Parameter(i).HasLimits()) {
               if (d > 0.5)
                  goto L26;
               d *= 10.;
               if (d > 0.5)
                  d = 0.51;
               continue;
            }
            d *= 10.;
         }

      L26:
         // get parameter name for i
         // (need separate scope for avoiding compl error when declaring name)
         print.Warn("2nd derivative zero for parameter", trafo.Name(trafo.ExtOfInt(i)),
                    "; MnHesse fails and will return diagonal matrix");

         for (unsigned int j = 0; j < n; j++) {
            double tmp = g2(j) < prec.Eps2() ? 1. : 1. / g2(j);
            vhmat(j, j) = tmp < prec.Eps2() ? 1. : tmp;
         }

         return MinimumState(st.Parameters(), MinimumError(vhmat, MinimumError::MnHesseFailed), st.Gradient(), st.Edm(),
                             mfcn.NumOfCalls());

      L30:
         double g2bfor = g2(i);
         g2(i) = 2. * sag / (d * d);
         grd(i) = (fs1 - fs2) / (2. * d);
         gst(i) = d;
         dirin(i) = d;
         yy(i) = fs1;
         double dlast = d;
         d = std::sqrt(2. * aimsag / std::fabs(g2(i)));
         if (trafo.Parameter(i).HasLimits())
            d = std::min(0.5, d);
         if (d < dmin)
            d = dmin;

         print.Debug("g1 =", grd(i), "g2 =", g2(i), "step =", gst(i), "d =", d, "diffd =", std::fabs(d - dlast) / d,
                     "diffg2 =", std::fabs(g2(i) - g2bfor) / g2(i));

         // see if converged
         if (std::fabs((d - dlast) / d) < strategy.HessianStepTolerance())
            break;
         if (std::fabs((g2(i) - g2bfor) / g2(i)) < strategy.HessianG2Tolerance())
            break;
         d = std::min(d, 10. * dlast);
         d = std::max(d, 0.1 * dlast);
      }
      vhmat(i, i) = g2(i);
      if (mfcn.NumOfCalls() > maxcalls) {

         // std::cout<<"maxcalls " << maxcalls << " " << mfcn.NumOfCalls() << "  " <<   st.NFcn() << std::endl;
         print.Warn("Maximum number of allowed function calls exhausted; will return diagonal matrix");

         for (unsigned int j = 0; j < n; j++) {
            double tmp = g2(j) < prec.Eps2() ? 1. : 1. / g2(j);
            vhmat(j, j) = tmp < prec.Eps2() ? 1. : tmp;
         }

         return MinimumState(st.Parameters(), MinimumError(vhmat, MinimumError::MnReachedCallLimit), st.Gradient(),
                             st.Edm(), mfcn.NumOfCalls());
      }
   }

   print.Debug("Second derivatives", g2);

   if (strategy.Strategy() > 0) {
      // refine first derivative
      HessianGradientCalculator hgc(mfcn, trafo, strategy);
      FunctionGradient gr = hgc(st.Parameters(), FunctionGradient(grd, g2, gst));
      // update gradient and step values
      grd = gr.Grad();
      gst = gr.Gstep();
   }

   // off-diagonal Elements
   // initial starting values
   if (n > 0) {
      unsigned int startParIndexOffDiagonal = 0;
      unsigned int endParIndexOffDiagonal = n * (n - 1) / 2;

      unsigned int offsetVect = 0;

      for (unsigned int in = startParIndexOffDiagonal; in < endParIndexOffDiagonal; in++) {

         int i = (in + offsetVect) / (n - 1);
         if ((in + offsetVect) % (n - 1) == 0) {
            offsetVect += i;
	 }
         int j = (in + offsetVect) % (n - 1) + 1;

         if ((i + 1) == j || in == startParIndexOffDiagonal) {
            x(i) += dirin(i);
         }

         x(j) += dirin(j);

         double fs1 = mfcn(x);
         double elem = (fs1 + amin - yy(i) - yy(j)) / (dirin(i) * dirin(j));
         vhmat(i, j) = elem;

         x(j) -= dirin(j);

         if (j % (n - 1) == 0 || in == endParIndexOffDiagonal - 1) {
            x(i) -= dirin(i);
	 }
      }
   }

   // verify if matrix pos-def (still 2nd derivative)

   print.Debug("Original error matrix", vhmat);

   MinimumError tmpErr = MnPosDef()(MinimumError(vhmat, 1.), prec);
   vhmat = tmpErr.InvHessian();

   print.Debug("PosDef error matrix", vhmat);

   int ifail = Invert(vhmat);
   if (ifail != 0) {

      print.Warn("Matrix inversion fails; will return diagonal matrix");

      MnAlgebraicSymMatrix tmpsym(vhmat.Nrow());
      for (unsigned int j = 0; j < n; j++) {
         double tmp = g2(j) < prec.Eps2() ? 1. : 1. / g2(j);
         tmpsym(j, j) = tmp < prec.Eps2() ? 1. : tmp;
      }

      return MinimumState(st.Parameters(), MinimumError(tmpsym, MinimumError::MnInvertFailed), st.Gradient(), st.Edm(),
                          mfcn.NumOfCalls());
   }

   FunctionGradient gr(grd, g2, gst);
   VariableMetricEDMEstimator estim;

   // if matrix is made pos def returns anyway edm
   if (tmpErr.IsMadePosDef()) {
      MinimumError err(vhmat, MinimumError::MnMadePosDef);
      double edm = estim.Estimate(gr, err);
      return MinimumState(st.Parameters(), err, gr, edm, mfcn.NumOfCalls());
   }

   // calculate edm for good errors
   MinimumError err(vhmat, 0.);
   double edm = estim.Estimate(gr, err);

   print.Debug("Hessian is ACCURATE. New state:", "\n  First derivative:", grd, "\n  Second derivative:", g2,
               "\n  Gradient step:", gst, "\n  Covariance matrix:", vhmat, "\n  Edm:", edm);

   return MinimumState(st.Parameters(), err, gr, edm, mfcn.NumOfCalls());
}

} // namespace TestStatistics
} // namespace RooFit

