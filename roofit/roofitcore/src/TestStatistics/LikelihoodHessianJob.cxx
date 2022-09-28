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

#include "ROOT/span.hxx"

#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/Numerical2PGradientCalculator.h"
#include "Minuit2/HessianGradientCalculator.h"
#include "Minuit2/MnPosDef.h"
#include "Minuit2/VariableMetricEDMEstimator.h"
#include "Minuit2/FCNBase.h"

#include "LikelihoodHessianJob.h"

namespace RooFit {
namespace TestStatistics {

using ROOT::Minuit2::FCNBase;
using ROOT::Minuit2::FunctionGradient;
using ROOT::Minuit2::HessianGradientCalculator;
using ROOT::Minuit2::MinimumError;
using ROOT::Minuit2::MinimumState;
using ROOT::Minuit2::MnAlgebraicSymMatrix;
using ROOT::Minuit2::MnAlgebraicVector;
using ROOT::Minuit2::MnFcn;
using ROOT::Minuit2::MnMachinePrecision;
using ROOT::Minuit2::MnPosDef;
using ROOT::Minuit2::MnPrint;
using ROOT::Minuit2::MnStrategy;
using ROOT::Minuit2::MnUserTransformation;
using ROOT::Minuit2::Numerical2PGradientCalculator;
using ROOT::Minuit2::VariableMetricEDMEstimator;

/// for implicit conversion from a raw double pointer to a matrix, only meant to add a simple double index operator
/// \note Very unsafe, no bounds checking!
/// operator() borrowed from LASymMatrix
struct RawDoublePtrMatrix {
   RawDoublePtrMatrix(double *data) : fData(data) {}
   double &operator()(unsigned int row, unsigned int col)
   {
      if (row > col)
         return fData[col + row * (row + 1) / 2];
      else
         return fData[row + col * (col + 1) / 2];
   }

private:
   double *fData;
};

/// \note trashsym will have its diagonal overwritten!
MinimumState diagonal_matrix_on_failure(enum MinimumError::Status failure_status, const MnFcn &mfcn,
                                        const MinimumState &st, const MnMachinePrecision &prec,
                                        const MnAlgebraicVector &g2, MnAlgebraicSymMatrix &trashsym)
{
   for (unsigned int j = 0; j < trashsym.Nrow(); j++) {
      double tmp = g2[j] < prec.Eps2() ? 1. : 1. / g2[j];
      trashsym(j, j) = tmp < prec.Eps2() ? 1. : tmp;
   }

   return MinimumState(st.Parameters(), MinimumError(trashsym, failure_status), st.Gradient(), st.Edm(),
                       mfcn.NumOfCalls());
}

void calculate_hessian_diagonal(const MnStrategy &strategy, const MnUserTransformation &trafo, unsigned int maxcalls,
                                double precision_eps2, double amin, double aimsag, const FCNBase *fcn,
                                RawDoublePtrMatrix vhmat, std::span<double> g2, std::span<double> gst,
                                std::span<double> grd, std::span<double> dirin, std::span<double> yy,
                                std::span<double> x);
void calculate_hessian_off_diagonal(double amin, std::span<double> dirin, std::span<double> yy, std::span<double> x,
                                    const FCNBase *fcn, RawDoublePtrMatrix vhmat);
/// Calculates the hessian from Minuit2 using LikelihoodHessianJob
///
/// This function can be used to replace MnHesse::calculator. Minuit2 will
/// then use it internally during minimization and during explicit calls to
/// calculate the hessian, e.g. from RooMinimizer::hesse().
MinimumState hessian_calculator(const MnStrategy &strategy, const MnFcn &mfcn, const MinimumState &st,
                                const MnUserTransformation &trafo, unsigned int maxcalls)
{
   MnPrint print("RooFit::MultiProcess enabled MnHesse calculator replacement");

   hessian_failure = MinimumError::Status::MnDefaultStatus;

   const MnMachinePrecision &prec = trafo.Precision();
   // make sure starting at the right place
   double amin = mfcn(st.Vec());
   double aimsag = std::sqrt(prec.Eps2()) * (std::fabs(amin) + mfcn.Up());

   // diagonal Elements first

   std::size_t n = st.Parameters().Vec().size();
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

   auto max_remaining_calls = maxcalls - mfcn.NumOfCalls();
   calculate_hessian_diagonal(strategy, trafo, max_remaining_calls, prec.Eps2(), amin, aimsag, &(mfcn.Fcn()),
                              vhmat.Data(), g2, gst, grd, dirin, yy, x);
   if (hessian_failure != MinimumError::Status::MnDefaultStatus) {
      return diagonal_matrix_on_failure(hessian_failure, mfcn, st, prec, g2, vhmat);
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
   calculate_hessian_off_diagonal(amin, dirin, yy, x, &(mfcn.Fcn()), vhmat.Data());

   // verify if matrix pos-def (still 2nd derivative)

   print.Debug("Original error matrix", vhmat);

   MinimumError tmpErr = MnPosDef()(MinimumError(vhmat, 1.), prec);
   vhmat = tmpErr.InvHessian();

   print.Debug("PosDef error matrix", vhmat);

   int ifail = Invert(vhmat);
   if (ifail != 0) {

      print.Warn("Matrix inversion fails; will return diagonal matrix");

      MnAlgebraicSymMatrix tmpsym(vhmat.Nrow());
      return diagonal_matrix_on_failure(MinimumError::MnInvertFailed, mfcn, st, prec, g2, tmpsym);
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
void calculate_hessian_off_diagonal(double amin, std::span<double> dirin, std::span<double> yy, std::span<double> x,
                                    const FCNBase *fcn, RawDoublePtrMatrix vhmat)
{ // off-diagonal Elements
   // initial starting values
   auto n = dirin.size();
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
            x[i] += dirin[i];
         }

         x[j] += dirin[j];

         double fs1 = (*fcn)(x);
         double elem = (fs1 + amin - yy[i] - yy[j]) / (dirin[i] * dirin[j]);
         vhmat(i, j) = elem;

         x[j] -= dirin[j];

         if (j % (n - 1) == 0 || in == endParIndexOffDiagonal - 1) {
            x[i] -= dirin[i];
         }
      }
   }
}

/// Return from this function when either the diagonal has been filled completely,
/// or when a failure has occurred. The hessian_failure flag can be set to the
/// appropriate failure condition. In hessian_calculator, this flag will trigger a
/// call to diagonal_matrix_on_failure to return the failure diagonal matrix.
///
/// \param strategy
/// \param trafo
/// \param maxcalls The maximum number of calls to fcn this function is allowed to make.
/// \param precision_eps2
/// \param amin
/// \param aimsag
/// \param mfcn
/// \param vhmat array with data for symmetric matrix of size n x n, so size in memory must be n (n + 1) / 2
/// \param g2 array of size n
/// \param gst array of size n
/// \param grd array of size n
/// \param dirin array of size n
/// \param yy array of size n
/// \param x array of size n
void calculate_hessian_diagonal(const MnStrategy &strategy, const MnUserTransformation &trafo, unsigned int maxcalls,
                                double precision_eps2, double amin, double aimsag, const FCNBase *fcn,
                                RawDoublePtrMatrix vhmat, std::span<double> g2, std::span<double> gst,
                                std::span<double> grd, std::span<double> dirin, std::span<double> yy,
                                std::span<double> &x)
{
   MnPrint print("RooFit::MultiProcess enabled MnHesse calculator replacement");

   std::size_t n_fcn_calls;

   auto n = g2.size();
   for (unsigned int i = 0; i < n; i++) {

      double xtf = x[i];
      double dmin = 8. * precision_eps2 * (fabs(xtf) + precision_eps2);
      double d = fabs(gst[i]);
      if (d < dmin)
         d = dmin;

      print.Debug("Derivative parameter", i, "d =", d, "dmin =", dmin);

      for (unsigned int icyc = 0; icyc < strategy.HessianNCycles(); icyc++) {
         double sag = 0.;
         double fs1 = 0.;
         double fs2 = 0.;
         for (unsigned int multpy = 0; multpy < 5; multpy++) {
            x[i] = xtf + d;
            fs1 = fcn(x);
            n_fcn_calls++;
            x[i] = xtf - d;
            fs2 = fcn(x);
            n_fcn_calls++;
            x[i] = xtf;
            sag = 0.5 * (fs1 + fs2 - 2. * amin);

            print.Debug("cycle", icyc, "mul", multpy, "\tsag =", sag, "d =", d);

            //  Now as F77 Minuit - check that sag is not zero
            if (sag != 0) {
               break;
            }
            if (trafo.Parameter(i).HasLimits()) {
               if (d > 0.5) {
                  print.Warn("2nd derivative zero for parameter", trafo.Name(trafo.ExtOfInt(i)),
                             "; MnHesse fails and will return diagonal matrix");

                  hessian_failure = MinimumError::MnHesseFailed;
                  return;
               }
               d *= 10.;
               if (d > 0.5)
                  d = 0.51;
               continue;
            }
            d *= 10.;
         }

         double g2bfor = g2[i];
         g2[i] = 2. * sag / (d * d);
         grd[i] = (fs1 - fs2) / (2. * d);
         gst[i] = d;
         dirin[i] = d;
         yy[i] = fs1;
         double dlast = d;
         d = sqrt(2. * aimsag / fabs(g2[i]));
         if (trafo.Parameter(i).HasLimits())
            d = std::min(0.5, d);
         if (d < dmin)
            d = dmin;

         print.Debug("g1 =", grd[i], "g2 =", g2[i], "step =", gst[i], "d =", d, "diffd =", fabs(d - dlast) / d,
                     "diffg2 =", fabs(g2[i] - g2bfor) / g2[i]);

         // see if converged
         if (fabs((d - dlast) / d) < strategy.HessianStepTolerance())
            break;
         if (fabs((g2[i] - g2bfor) / g2[i]) < strategy.HessianG2Tolerance())
            break;
         d = std::min(d, 10. * dlast);
         d = std::max(d, 0.1 * dlast);
      }
      vhmat(i, i) = g2[i];
      if (n_fcn_calls > maxcalls) {
         print.Warn("Maximum number of allowed function calls exhausted; will return diagonal matrix");

         hessian_failure = MinimumError::MnReachedCallLimit;
         return;
      }
   }
}

} // namespace TestStatistics
} // namespace RooFit
