// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann, E.G.P. Bos    2013-2018
//
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 *                                                                    *
 **********************************************************************/
/*
 * NumericalDerivator.h
 *
 *  Original version created on: Aug 14, 2013
 *      Authors: L. Moneta, J. T. Offermann
 *  Modified version created on: Sep 27, 2017
 *      Author: E. G. P. Bos
 */

#ifndef ROOT_Minuit2_NumericalDerivator
#define ROOT_Minuit2_NumericalDerivator

#ifndef ROOT_Math_IFunctionfwd
#include <Math/IFunctionfwd.h>
#endif

#include <vector>
#include "Fit/ParameterSettings.h"
#include "Minuit2/SinParameterTransformation.h"
#include "Minuit2/SqrtUpParameterTransformation.h"
#include "Minuit2/SqrtLowParameterTransformation.h"
#include "Minuit2/MnMachinePrecision.h"

namespace ROOT {
namespace Minuit2 {

// Holds all necessary derivatives and associated numbers (per parameter) used in the NumericalDerivator class.
struct DerivatorElement {
   double derivative;
   double second_derivative;
   double step_size;
};

class NumericalDerivator {
public:
   explicit NumericalDerivator(bool always_exactly_mimic_minuit2 = true);
   NumericalDerivator(const NumericalDerivator &other);
   NumericalDerivator(double step_tolerance, double grad_tolerance, unsigned int ncycles, double error_level,
                             bool always_exactly_mimic_minuit2 = true);
   virtual ~NumericalDerivator();

   void setup_differentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *cx,
                            const std::vector<ROOT::Fit::ParameterSettings> &parameters);
   std::vector<DerivatorElement>
   Differentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                 const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                 const std::vector<DerivatorElement>& previous_gradient);

   DerivatorElement
   partial_derivative(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                      const std::vector<ROOT::Fit::ParameterSettings> &parameters, unsigned int i_component, DerivatorElement previous);
   DerivatorElement fast_partial_derivative(const ROOT::Math::IBaseFunctionMultiDim *function,
                                const std::vector<ROOT::Fit::ParameterSettings> &parameters, unsigned int i_component, const DerivatorElement& previous);
   DerivatorElement operator()(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                                                 const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                                 unsigned int i_component, const DerivatorElement& previous);

   double GetFValue() const { return fVal; }
   void SetStepTolerance(double value);
   void SetGradTolerance(double value);
   void SetNCycles(int value);

   double Int2ext(const ROOT::Fit::ParameterSettings &parameter, double val) const;
   double Ext2int(const ROOT::Fit::ParameterSettings &parameter, double val) const;
   double DInt2Ext(const ROOT::Fit::ParameterSettings &parameter, double val) const;

   void SetInitialGradient(const ROOT::Math::IBaseFunctionMultiDim *function,
                           const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                           std::vector<DerivatorElement>& gradient);

   void set_step_tolerance(double step_tolerance);
   void set_grad_tolerance(double grad_tolerance);
   void set_ncycles(unsigned int ncycles);
   void set_error_level(double error_level);

private:
   double fStepTolerance = 0.5;
   double fGradTolerance = 0.1;
   unsigned int fNCycles = 2;
   double Up = 1;
   double fVal = 0;

   std::vector<double> vx, vx_external;
   double dfmin;
   double vrysml;

   // MODIFIED: Minuit2 determines machine precision in a slightly different way than
   // std::numeric_limits<double>::epsilon()). We go with the Minuit2 one.
   ROOT::Minuit2::MnMachinePrecision precision;

   ROOT::Minuit2::SinParameterTransformation fDoubleLimTrafo;
   ROOT::Minuit2::SqrtUpParameterTransformation fUpperLimTrafo;
   ROOT::Minuit2::SqrtLowParameterTransformation fLowerLimTrafo;

private:
   bool _always_exactly_mimic_minuit2;

public:
   bool always_exactly_mimic_minuit2() const;
   void set_always_exactly_mimic_minuit2(bool flag);

private:
   std::vector<double> vx_fVal_cache;
};

std::ostream& operator<<(std::ostream& out, const DerivatorElement value);

} // namespace Minuit2
} // namespace ROOT

#endif  // ROOT_Minuit2_NumericalDerivator