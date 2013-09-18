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


#include <math.h>
#include <stdio.h>

#include "SpecialFunctionGenerator.h"

SpecialFunctionGenerator::SpecialFunctionGenerator()
{

#if (DEBUG_SF == 1)
  branch1 = 0;
  branch2 = 0;
#endif
}


SpecialFunctionGenerator::~SpecialFunctionGenerator()
{
#if (DEBUG_SF == 1)
  fprintf(stderr,"branch 1: %ld branch 2: %ld\n",branch1,branch2);
#endif
}


double SpecialFunctionGenerator::BesselJIntegralOrder(int order, double x)
{
   double result=0;
   if (order == 0) {
      result = BesselJ0(x);
   }
   else if (order == 1) {
      result = BesselJ1(x);
   }
   else if (order > 1) {

   }
   else if (order < 0){

   }
  
   return result;
}

// this is similar to that found in Numerical Recipes
double SpecialFunctionGenerator::BesselJ0(double x)
{
   double ax = fabs(x);
   double result = 0;

   if (ax < 8.0) {
      double r1,r2;
      double y = x*x;
      r1 = 57568490571.0 + y*(-13362590354.0 + y*
			      (651619640.7 + y*
			       (-11214424.18+y*
				(77392.33017+y*
				 (-184.9052456)))));
	
      r2 = 57568490411.0 + y*(1029532685.0 + y*
			      (9494680.718 +y*
			       (59272.64853+y*
				(267.8532712+y))));

      result = r1 / r2; 
   } 
   else {
      double z = 8.0 / ax;
      double y = z * z;
      double xx = ax - 0.785398164;
      double r1,r2;

      r1 = 1.0 + y*(-0.1098628627e-2 + y*
		    (0.2734510407e-4 + y*
		     (-0.2073370639e-5 + y*
		      0.2093887211e-6)));

      r2 = -0.1562499995e-1 + y*
	 (0.1430488765e-3 + y*
	  (-0.6911147651e-5 + y*
	   (0.7621095161e-6 - y*
	    0.93494512e-7)));

      result = sqrt(0.636619772/ax) * (cos(xx)*r1-z*sin(xx)*r2);
      
   }

   return result; 
}


// this is similar to that found in Numerical Recipes
double SpecialFunctionGenerator::BesselJ1(double x)
{
   double ax = fabs(x);
   double result = 0;

   if (ax < 8.0) {
#if (DEBUG_SF == 1)
     branch1++;
#endif
      double r1,r2;
      double y = x*x;
      r1 = x*(72362614232.0 + y*(-7895059235.0 + y*
				 (242396853.1 + y*
				  (-2972611.439+y*
				   (15704.48260+y*
				    (-30.16036606))))));
	
      r2 = 144725228442.0 + y*(2300535178.0 + y*
			      (18583304.74 +y*
			       (99447.43394+y*
				(376.9991397+y))));

      result = r1 / r2; 
   } 
   else {
#if (DEBUG_SF == 1)
     branch2++;
#endif
      double z = 8.0 / ax;
      double y = z * z;
      double xx = ax - 2.356194491;
      double r1,r2;

      r1 = 1.0 + y*(0.183105e-2 + y*
		    (-0.3516396496e-4 + y*
		     (0.2457520174e-5 + y*
		      (-0.240337019e-6))));

      r2 = 0.04687499995 + y*
	 (-0.2002690873e-3 + y*
	  (0.8449199096e-5 + y*
	   (-0.88228987e-6 - y*
	    0.105787412e-6)));

      result = sqrt(0.636619772/ax) * (cos(xx)*r1-z*sin(xx)*r2);
      if (x < 0) result = -result;
   }

   return result;    
}
