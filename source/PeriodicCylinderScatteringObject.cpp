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

#include "HS3DVector.h"
#include "HSRandomStream.h"
#include "PeriodicCylinderScatteringObject.h"
#include "SpecialFunctionGenerator.h"


vector<HSComplex> PeriodicCylinderScatteringObject::CalculateScatteredAmplitude(
   vector<double> &kxArray, vector<double> &kyArray,vector<double> &kzArray)
{   
   vector<HSComplex> disc_factor(kxArray.size(),0);
   vector<HSComplex> axial_factor(kxArray.size(),1);


   calculate_disc_structure_factor(disc_factor,kxArray,kyArray,kzArray);

   calculate_axial_structure_factor(axial_factor,kxArray,kyArray,kzArray);


   vector<HSComplex> result(kxArray.size());

   

   for (long i=0;i<result.size();i++) {
      double phase = x * kxArray[i] + y * kyArray[i] + z * kzArray[i];

      result[i] = disc_factor[i] * axial_factor[i] * HSComplex(cos(phase),sin(phase));
   }

   return result;
}


// Calculates the scattering by placing numpoints scattering centres inside the cylinder.
vector<HSComplex> PeriodicCylinderScatteringObject::CalculateScatteredAmplitudeRandomized(
    vector<double> &kxArray, vector<double> &kyArray,
    vector<double> &kzArray, long seed, long numpoints)
{
    HSRandomStream rstream(seed);


    vector<HSComplex> result(kzArray.size(),0);
    long i = 0;
    double R2 = R*R;
    while (i < numpoints) {
	double x = rstream.continuous_uniform(-1,1) * R;
	double y = rstream.continuous_uniform(-1,1) * R;
	double z = rstream.continuous_uniform(0,1) * L * n;

	if ((x*x+y*y) <= R2) {
	    i++;

	    double rho = 1;

	    if (((x-floor(x/L)*L) / L) < sigma) {
		rho = rho1;
	    }
	    else {
		rho = rho2;
	    }

	    double zp = -sin(alpha)*x+cos(alpha) * z;
	    double xp = x * cos(alpha) * cos(beta) + z * sin(alpha) * cos(beta) - y * sin(beta);
	    double yp = x * cos(alpha) * sin(beta) + z * sin(alpha) * sin(beta) + y * cos(beta);
	    

	    for (long j=0;j<kxArray.size();j++) {
		double phase = xp * kxArray[j] + yp * kyArray[j] + zp * kzArray[j];

		result[j] += HSComplex(rho*cos(phase),rho*sin(phase));
	    }

	}

    }


    return result;
}




ScatteringObject* PeriodicCylinderScatteringObject::InstantiateObject(char *s)
{
   vector<double> params;

   double temp;
   long n=0;
   long n2;
   while (sscanf(s+n,"%lf %n",&temp,&n2) == 1) {
      params.push_back(temp);
      n = n + n2;
   }

   if (params.size() != 11) {
      // some problems with the string
      fprintf(stderr,"problems with string %s\n",s);
      return NULL;
   }
   else {
      PeriodicCylinderScatteringObject *obj = new PeriodicCylinderScatteringObject();

      obj->x = params[0];
      obj->y = params[1];
      obj->z = params[2];
      obj->R = params[3];
      obj->alpha = params[4] / 180.0 * M_PI; // convert angle to radians
      obj->beta = params[5] / 180.0 * M_PI; // convert angle to radians
      obj->L = params[6];
      obj->n = (long)params[7];
      obj->rho1 = params[8];
      obj->rho2 = params[9];
      obj->sigma = params[10];

      return obj;
   }
}


void PeriodicCylinderScatteringObject::calculate_disc_structure_factor(
   vector<HSComplex> &outStructureFactor,vector<double> &kxArray, 
   vector<double> &kyArray,vector<double> &kzArray)
{
   double c_alpha = cos(alpha);
   double c_beta = cos(beta);
   double s_alpha = sin(alpha);
   double s_beta = sin(beta);


   SpecialFunctionGenerator gen;

   for (long i=0;i<kxArray.size();i++) {
      double kx = kxArray[i];
      double ky = kyArray[i];
      double kz = kzArray[i];
 
      double kxp = c_alpha * (kx*c_beta + ky*s_beta) - kz*s_alpha; 
      double kyp = ky*c_beta - kx*s_beta;
     
      double a = sqrt(kxp*kxp+kyp*kyp);

      if (a != 0) {
	 double result = 2*M_PI*R/a*gen.BesselJIntegralOrder(1,R*a);
	 outStructureFactor[i].set(result,0);
      }
      else {
	 outStructureFactor[i].set(M_PI*R*R,0);
      }
   }
}






void PeriodicCylinderScatteringObject::calculate_axial_structure_factor(
   vector<HSComplex> &outStructureFactor,vector<double> &kxArray, 
   vector<double> &kyArray,vector<double> &kzArray)
{
   HS3DVector ba(sin(alpha)*cos(beta),sin(alpha)*sin(beta),cos(alpha));

   // precalculate the trigonometric functions
   double c_alpha = cos(alpha);
   double s_alpha = sin(alpha);
   double c_beta = cos(beta);
   double s_beta = sin(beta);

   double re=0,im=0;

   for (long i=0;i<kxArray.size();i++) {
      double kx = kxArray[i];
      double ky = kyArray[i];
      double kz = kzArray[i];

      double kzp = s_alpha * (kx*c_beta + ky*s_beta) + kz*c_alpha;



      double sin_Lk = sin(L*kzp);
      double cos_Lk = cos(L*kzp);
      double cos_n1Lk = cos((n-1)*L*kzp);
      double sin_n1Lk = sin((n-1)*L*kzp);
      double cos_nLk = cos_n1Lk*cos_Lk - sin_n1Lk*sin_Lk;
      double sin_nLk = cos_n1Lk*sin_Lk + cos_Lk*sin_n1Lk;
      double sin_Lk_sigma = sin(L*kzp*sigma);
      double cos_Lk_sigma = cos(L*kzp*sigma);

      if ((kzp != 0) && ((2-2*cos_Lk) != 0)) {
	 double c1 = 1 / (kzp * (2 - 2 * cos_Lk));
	 double re1 = c1 * (rho1 * sin_Lk_sigma + rho2 * sin_Lk - rho2 * sin_Lk_sigma);
	 double im1 = -c1 * (rho1 * cos_Lk_sigma - rho1 + rho2 * cos_Lk - rho2 * cos_Lk_sigma);

	 double re2 = (1 - cos_Lk - cos_nLk + cos_n1Lk);
	 double im2 = -(-sin_Lk + sin_nLk - sin_n1Lk);



	 re = re1 * re2 - im1 * im2;
	 im = re1 * im2 + im1 * re2;
      }
      else { 
	 re = n * L * (rho1 * sigma + (1 - sigma) * rho2);
	 im = 0;
      }

      outStructureFactor[i].set(re,im);
   }
}
