// @(#)root/mathmore:$Id$
// Author: L. Moneta June 2009

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/


// Header file for class MinimTransformFunction

#ifndef ROOT_Math_MinimTransformFunction
#define ROOT_Math_MinimTransformFunction


#include "Math/IFunction.h"

#include "Math/MinimTransformVariable.h"


#include <vector>
#include <map>

namespace ROOT {

   namespace Math {



/**
   MinimTransformFunction class to perform a transformations on the
   variables to deal with fixed or limited variables (support both double and single bounds)
   The class manages the passed function pointer

   @ingroup MultiMin
*/
class MinimTransformFunction : public IMultiGradFunction {

public:

   typedef  ROOT::Math::IMultiGradFunction BaseGradFunc;
   typedef  ROOT::Math::IMultiGradFunction::BaseFunc BaseFunc;


   /**
     Constructor from a IMultiGradFunction interface (which is managed by the class)
     vector specifying the variable types (free, bounded or fixed, defined in enum EMinimVariableTypes )
     variable values (used for the fixed ones) and a map with the bounds (for the bounded variables)

   */
   MinimTransformFunction ( const IMultiGradFunction * f, const std::vector<ROOT::Math::EMinimVariableType> & types, const std::vector<double> & values,
                            const std::map<unsigned int, std::pair<double, double> > & bounds);


   /**
      Destructor (delete function pointer)
   */
   ~MinimTransformFunction ()  {
      if (fFunc) delete fFunc;
   }


   // method inherited from IFunction interface

   unsigned int NDim() const { return fIndex.size(); }

   unsigned int NTot() const { return fFunc->NDim(); }

   /// clone:  not supported (since unique_ptr used in the fVariables)
   IMultiGenFunction * Clone() const {
      return 0;
   }


   /// transform from internal to external
   /// result is cached also inside the class
   const std::vector<double> & Transformation( const std::vector<double> & x) const {
      Transformation(x, fX);
      return fX;
  }


   /// transform from internal to external
   void Transformation( const std::vector<double> & xint, std::vector<double> & xext) const;

   /// inverse transformation (external -> internal)
   void  InvTransformation(const std::vector<double> & xext,  std::vector<double> & xint) const;

   /// inverse transformation for steps (external -> internal) at external point x
   void  InvStepTransformation(const std::vector<double> & x, const std::vector<double> & sext,  std::vector<double> & sint) const;

   ///transform gradient vector (external -> internal) at internal point x
   void GradientTransformation(const std::vector<double> & x, const std::vector<double> &gExt, std::vector<double> & gInt) const;

   ///transform covariance matrix (internal -> external) at internal point x
   /// use row storages for matrices  m(i,j) = rep[ i * dim + j]
   void MatrixTransformation(const std::vector<double> & x, const std::vector<double> &covInt, std::vector<double> & covExt) const;

   // return original function
   const IMultiGradFunction *OriginalFunction() const { return fFunc; }


private:

   /// function evaluation
   virtual double DoEval(const std::vector<double> & x) const {
#ifndef DO_THREADSAFE
      return (*fFunc)(Transformation(x));
#else
      std::vector<double> xext(fVariables.size() );
      Transformation(x, &xext[0]);
      return (*fFunc)(&xext[0]);
#endif
   }

   /// calculate derivatives
   virtual double DoDerivative (const std::vector<double> & x, unsigned int icoord  ) const {
      const MinimTransformVariable & var = fVariables[ fIndex[icoord] ];
      double dExtdInt = (var.IsLimited() ) ? var.DerivativeIntToExt( x[icoord] ) : 1.0;
      double deriv =  fFunc->Derivative( Transformation(x) , fIndex[icoord] );
      //std::cout << "Derivative icoord (ext)" << fIndex[icoord] << "   dtrafo " << dExtdInt << "  " << deriv << std::endl;
      return deriv * dExtdInt;
   }

   // copy constructor for this class (disable by having it private)
   MinimTransformFunction( const MinimTransformFunction & ) :
      BaseFunc(), BaseGradFunc()
   {}

   // assignment operator for this class (disable by having it private)
   MinimTransformFunction & operator= ( const MinimTransformFunction & ) {
      return *this;
   }



private:

   mutable std::vector<double>  fX;              // internal cached of external values
   std::vector<MinimTransformVariable> fVariables;    // vector of variable settings and tranformation function
   std::vector<unsigned int>      fIndex;        // vector with external indices for internal variables
   const IMultiGradFunction * fFunc;             // user function

};

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_MinimTransformFunction */
