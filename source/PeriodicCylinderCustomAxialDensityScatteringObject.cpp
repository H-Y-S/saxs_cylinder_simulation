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
#include "PeriodicCylinderCustomAxialDensityScatteringObject.h"
#include "SpecialFunctionGenerator.h"

vector<HSComplex> PeriodicCylinderCustomAxialDensityScatteringObject::axial_scattering_amplitude;
double PeriodicCylinderCustomAxialDensityScatteringObject::k_step;
vector<HSComplex> PeriodicCylinderCustomAxialDensityScatteringObject::CalculateScatteredAmplitude(
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


// Calculates the scattering by placing numpoints scattering centres
// inside the cylinder.
vector<HSComplex> PeriodicCylinderCustomAxialDensityScatteringObject::CalculateScatteredAmplitudeRandomized(
    vector<double> &kxArray, vector<double> &kyArray,
    vector<double> &kzArray, long seed, long numpoints)
{
  fprintf(stderr,"PeriodicCylinderCustomAxialDensityScatteringObject does not support randomized calculation\n");
  vector<HSComplex> result(1);
   return result;;
}




ScatteringObject* PeriodicCylinderCustomAxialDensityScatteringObject::InstantiateObject(char *s)
{
  if (axial_scattering_amplitude.size() == 0) {
    fprintf(stderr,"Cannot instantiate object due to lacking axial electron density file\n");
    exit(1);
  }

  vector<double> params;

  double temp;
  long n=0;
  long n2;
  while (sscanf(s+n,"%lf %n",&temp,&n2) == 1) {
    params.push_back(temp);
    n = n + n2;
  }

  if (params.size() != 8) {
    // some problems with the string
    fprintf(stderr,"problems with string %s\n",s);
    return NULL;
  }
  else {
    PeriodicCylinderCustomAxialDensityScatteringObject *obj = 
      new PeriodicCylinderCustomAxialDensityScatteringObject();

    obj->x = params[0];
    obj->y = params[1];
    obj->z = params[2];
    obj->R = params[3];
    obj->alpha = params[4] / 180.0 * M_PI; // convert angle to radians
    obj->beta = params[5] / 180.0 * M_PI; // convert angle to radians
    obj->L = params[6];
    obj->n = (long)params[7];
    
    return obj;
  }
}


void PeriodicCylinderCustomAxialDensityScatteringObject::SetAxialDensity(const char *filename)
{
  FILE *fp = fopen(filename,"r");


  if (fp != NULL) {
    vector<double> coord;
    vector<double> density;
    
    // first read the file
    float c,d;
    while (fscanf(fp,"%f %f",&c,&d) == 2) {
      coord.push_back(c);
      density.push_back(d);
    }


    // calculate the scattered amplitude at sufficient density and to
    // large enough values of k*L (let us cut this at 1000 for these objects 
    // for objects of 1000 Å in length at k = 1/Å) because after that things 
    // are not really SAXS anyway.
    InitializeAxialScatteringAmplitude(coord,density);
  }
  else {
    fprintf(stderr,"Axial density file: %s could not be opened\n",filename);
    exit(1);
  }
}



void PeriodicCylinderCustomAxialDensityScatteringObject::calculate_disc_structure_factor(
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






void PeriodicCylinderCustomAxialDensityScatteringObject::calculate_axial_structure_factor(
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


      // calculate the projection length of k-vector to the cylinder axis.
      double kzp = s_alpha * (kx*c_beta + ky*s_beta) + kz*c_alpha;
      
      double kzp_normalized = fabs(kzp) * L / MAX_PERIOD_LENGTH;
      kzp = fabs(kzp);
    
      // calculate the k-vector index
      double index = kzp_normalized / k_step;
      long low_ind = (long)floor(index);
      long high_ind = (long)ceil(index);
  
      HSComplex low_val = axial_scattering_amplitude[low_ind];
      HSComplex high_val = axial_scattering_amplitude[high_ind];

      double re1;
      double im1;

      if (low_ind != high_ind) {
	re1 = low_val.Re() * (high_ind-index) + 
	  high_val.Re() * (index-low_ind);

	im1 = low_val.Im() * (high_ind-index) + 
	  high_val.Im() * (index-low_ind);
      }
      else {
	re1 = low_val.Re();
	im1 = low_val.Im();
      }

      double sin_Lk = sin(L*kzp);
      double cos_Lk = cos(L*kzp);
      double cos_n1Lk = cos((n-1)*L*kzp);
      double sin_n1Lk = sin((n-1)*L*kzp);
      double cos_nLk = cos_n1Lk*cos_Lk - sin_n1Lk*sin_Lk;
      double sin_nLk = cos_n1Lk*sin_Lk + cos_Lk*sin_n1Lk;

      double re2,im2,re,im;

      if ((2-2*cos_Lk) != 0) {
	double c1 = 1 / (2 - 2 * cos_Lk);
	re2 = c1 * (1 - cos_Lk - cos_nLk + cos_n1Lk);
	im2 = -c1 * (-sin_Lk + sin_nLk - sin_n1Lk);
	
	re = re1 * re2 - im1 * im2;
	im = re1 * im2 + im1 * re2;
      }
      else {
	re = axial_scattering_amplitude[0].Re()*n;
	im = 0;
      }

      outStructureFactor[i].set(re,im);
   }
}





void PeriodicCylinderCustomAxialDensityScatteringObject::InitializeAxialScatteringAmplitude
(const vector<double> &coord, const vector<double> &dens)
{
  // the smallest required k-step is inversely proportional to the
  // period length
  float max_c = coord[coord.size()-1];
  k_step = 1 / (10.0*MAX_PERIOD_LENGTH);
  
  for (float k=0;k<=1;k += k_step) {
    HSComplex v = 0;

    // coords are in relative units from 0 to 1. 
    for (long i=0;i<coord.size();i++) {
      float c = coord[i] * MAX_PERIOD_LENGTH;
      float d = dens[i];

      HSComplex b(d * cos(k*c),d * sin(k*c));
      v = v + b;
    }

    axial_scattering_amplitude.push_back(v);
  }
}
											   
