/*
  Author: Heikki Suhonen (heikki.suhonen@helsinki.fi)
  
  This source file is a part of a SAXS simulation package made at the
  University of Helsinki, Department of Physical Sciences, Division of
  X-ray Physics. The permission to use this package for academic
  research purposes is granted. If you publish any results that
  include results made with this package we would like you to include
  a citation to H. Suhonen et al. Phys. Med. Biol. 50 (2005),
  5401-5416. For other possible uses of the file, please contact the
  author.

  The source file is provided "as is" and its correctness is not
  guaranteed. If you find any suspicious parts in the code, please
  inform the author.
*/


#ifndef _SPECIAL_FUNCTION_GENERATOR_H
#define _SPECIAL_FUNCTION_GENERATOR_H


#define DEBUG_SF 0

class SpecialFunctionGenerator {
  public:
   SpecialFunctionGenerator();
   ~SpecialFunctionGenerator();

   double BesselJIntegralOrder(int order, double x);
  

   inline   double HermitePolynomial(long order, double x);



  private:
   double BesselJ0(double x);
   double BesselJ1(double x);

#if (DEBUG_SF == 1)
   long branch1;
   long branch2;
#endif
};


// Recursive Hermite polynomial, very slow for large order, seems to work OK 
// at least for the first four orders.
inline double SpecialFunctionGenerator::HermitePolynomial(long order, double x)
{
   if (order == 0) {
      return 1;
   }
   else if (order == 1) {
      return 2*x;
   }
   else if (order > 1) {
      return 2*x*HermitePolynomial(order-1,x) - 2*(order-1)*HermitePolynomial(order-2,x);
   }
}



#endif // _SPECIAL_FUNCTION_GENERATOR_H
