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

#include <algorithm>
#include <curses.h>
#include <math.h>
#include <stdio.h>
#include <vector>


#include "HS3DVector.h"
#include "HSRandomStream.h"
#include "Utilities.h"

using namespace std;


vector<HS3DVector> unit_lattice;

void rotate_lattice(vector<HS3DVector> &lattice, double angle);
void tilt_lattice(vector<HS3DVector> &lattice, double alpha, double beta);
void read_param(char *filename);

struct interval {
      double begin;
      double end;

      bool operator<(const interval &b) const {
	 return begin < b.begin;
      }
};



long unit_count = 1;
long cylinders_per_unit_min = 1;
long cylinders_per_unit_max = 30;

double R_avg = 300;
double R_stdev = 20;
double L_avg = 670;
double L_stdev = 0;
double zdisp = 1000; // displacement within single bunch

double alpha_stdev = 0; 

double rho1 = 0;
double rho2 = 1.0;
double period_proportion = 0.5;
double period_proportion_stdev = 0;
long period_min = 10;
long period_max = 10;

bool do_tilt_lattice = false;
bool use_random_orientation = false;

double packdist_radius_fraction = 1.0;
double packdist_addition = 0;

char param_string[8192]; 


int main(int argc, char *argv[])
{
   if (argc > 1) {
      read_param(argv[1]);
   }
   else {
      fprintf(stderr,"Warning, using default parameters\n");
   }

   HSRandomStream rand(1),rand2(2),rand3(3),rand4(5),rand5(7),rand6(11),rand7(13);

   vector<HS3DVector> final_lattice;
   vector<double> final_radii;
   vector<double> final_lengths;


   printf("###\n");
   printf("%s",param_string);
   printf("###\n");


   for (long i=0;i<unit_count;i++) {
      long cylinder_count = rand.discrete_uniform(cylinders_per_unit_min,cylinders_per_unit_max);

      final_lattice.clear();
      final_radii.clear();
      final_lengths.clear();



      double alpha = rand3.gaussian_distribution(0,alpha_stdev);
      double beta = rand3.continuous_uniform(0,360);

      vector<HS3DVector> lattice;

      if (use_random_orientation) {
	 double t = rand2.continuous_uniform(0,1);
	 alpha = acos(t) / M_PI * 180;
      }


      double limit_rad = 0;

      for (long k=0;k<cylinder_count;k++) {
	 double R = rand2.gaussian_distribution(R_avg,R_stdev);

	 double alpha1 = rand3.continuous_uniform(0,2*M_PI);
	 double x = limit_rad * cos(alpha1);
	 double y = limit_rad * sin(alpha1);
	 
	 HS3DVector p1(x,y,0);

	 double mind = 1e10;
	 long minj = -1;


	 for (long j=0;j<final_lattice.size();j++) {
	    HS3DVector p2 = final_lattice[j];
	    
	    double dist = p1.dist(p2);

	    if (dist < mind) {
	       mind = dist;
	       minj = j;
	    }
	 }

	 if (minj >= 0) {
	    HS3DVector p2 = final_lattice[minj];
	    double R2 = final_radii[minj];
	    
	    double minbeta_neg = M_PI;
	    double minbeta_pos = M_PI;

	    HS3DVector a = p2;
	    HS3DVector a_hat = p2 * (1.0 / p2.length());
	    HS3DVector an_hat(a_hat.Y(),-a_hat.X(),0);

	    for (long j=0;j<final_lattice.size();j++) {
	       if (j != minj) {
		  HS3DVector c = final_lattice[j];
		  HS3DVector b = c - a;

		  double rand_add = rand7.continuous_uniform(0,packdist_addition);

		  double d1 = (R + R2 + rand_add) * packdist_radius_fraction;
		  double d2 = (R + final_radii[j]) * packdist_radius_fraction;
		  double d3 = b.length();

		  if (d3 < d1+d2) {
		     long sign = long(a.X()*c.Y()-a.Y()*c.X());
		     
		     double arg = a.dot_product(b) / (a.length() * b.length());
		     if (fabs(arg) > 1.0) {
			arg = arg / fabs(arg);
		     }

		     double gamma = acos(arg);
		     double alfa = acos((d3*d3+d1*d1-d2*d2)/(2*d3*d1));

		     if (sign >= 0) {
			minbeta_pos = min_c(minbeta_pos,gamma-alfa);

			if (gamma+alfa > M_PI) {
			   minbeta_neg = min_c(minbeta_neg,2*M_PI-gamma-alfa);
			}
		     }
		     else {
			minbeta_neg = min_c(minbeta_neg,gamma-alfa);

			if (gamma+alfa > M_PI) {
			   minbeta_pos = min_c(minbeta_pos,2*M_PI-gamma-alfa);
			}
		     }
		     if (isnan(gamma) || isnan(alfa)) {
			fprintf(stderr,"isnan: %f %f, a*b = %f, |a| = %f, |b| = %f\n",gamma,alfa,a.dot_product(b),
				a.length(),b.length());
		     }
		  }
	       }
	    }


	    if (minbeta_neg > minbeta_pos) {
	       p1 = a + a_hat*(R+R2+packdist_addition)*packdist_radius_fraction*cos(minbeta_neg) + 
		  an_hat*(R+R2+packdist_addition)*packdist_radius_fraction*sin(minbeta_neg); 
	    }
	    else {
	       p1 = a + a_hat*(R+R2+packdist_addition)*packdist_radius_fraction*cos(minbeta_pos) - 
		  an_hat*(R+R2+packdist_addition)*packdist_radius_fraction*sin(minbeta_pos); 
	    }
	    final_lattice.push_back(p1);


	    limit_rad = max_c(limit_rad,p1.length()+R);
	 }
	 else {
	    final_lattice.push_back(HS3DVector(0,0.1,0));
	    limit_rad = R;
	 }

	 final_radii.push_back(R);
	 final_lengths.push_back(rand4.gaussian_distribution(L_avg,L_stdev));
      }


      // Add the axial displacement for fibrils
      for (long j=0;j<final_lattice.size();j++) {
	 HS3DVector p1 = final_lattice[j];
	 HS3DVector p2(p1.X(),p1.Y(),rand3.continuous_uniform(0,zdisp));

	 final_lattice[j] = p2;
      }
      




      if (do_tilt_lattice == 1)
	 tilt_lattice(final_lattice,alpha / 180 * M_PI,beta / 180 * M_PI);


      // the scattering unit is done with the values in final_... vectors
      for (long j=0;j<final_lattice.size();j++) {
	 double x = final_lattice[j].X();
	 double y = final_lattice[j].Y();
	
	 printf("cylper %f %f %f %f %f %f %f %ld %f %f %f\n",final_lattice[j].X(),
		final_lattice[j].Y(),final_lattice[j].Z(),final_radii[j],
		alpha,beta,final_lengths[j],
		rand5.discrete_uniform(period_min,period_max),
		rho1,rho2,
		period_proportion+rand5.gaussian_distribution(0,period_proportion_stdev));
      }

      if (i < unit_count-1)
	 printf("*\n");
   }


   return 0;
}


void rotate_lattice(vector<HS3DVector> &lattice, double angle)
{
   for (int i=0;i<lattice.size();i++) {
      double x,y,z;
      x = lattice[i].X() * cos(angle) - lattice[i].Y() * sin(angle);
      y = lattice[i].X() * sin(angle) + lattice[i].Y() * cos(angle);
      z = lattice[i].Z();

      lattice[i] = HS3DVector(x,y,z);
   }
}


void tilt_lattice(vector<HS3DVector> &lattice, double alpha, double beta)
{
   for (int i=0;i<lattice.size();i++) { 
      HS3DVector p = lattice[i];
      double x,y,z;
      
      x = (cos(alpha) * p.X() + sin(alpha) * p.Z()) * cos(beta) - sin(beta) * p.Y();
      y = cos(beta) * p.Y() + (cos(alpha) * p.X() + sin(alpha) * p.Z()) * sin(beta);
      z = cos(alpha) * p.Z() - sin(alpha) * p.X();


      lattice[i] = HS3DVector(x,y,z);
   }
}




void read_param(char *filename)
{
   FILE *fp = fopen(filename,"r");
   char option[1024];
   float value;
   if (fp != NULL) {
      while (fscanf(fp,"%s = %f\n",option,&value) > 0) {
	 if (strcmp(option,"unit_count") == 0) {
	    unit_count = (long)value;
	 }
	 else if (strcmp(option,"cylinders_per_unit_min") == 0) {
	    cylinders_per_unit_min = long(value);
	 }
	 else if (strcmp(option,"cylinders_per_unit_max") == 0) {
	    cylinders_per_unit_max = long(value);
	 }
	 else if (strcmp(option,"r_avg") == 0) {
	    R_avg = value;
	 }
	 else if (strcmp(option,"r_stdev") == 0) {
	    R_stdev = value;
	 }
	 else if (strcmp(option,"l_avg") == 0) {
	    L_avg = value;
	 }
	 else if (strcmp(option,"l_stdev") == 0) {
	    L_stdev = value;
	 }
	 else if (strcmp(option,"alpha_stdev") == 0) {
	    alpha_stdev = value;
	 }
	 else if (strcmp(option,"zdisp") == 0) {
	    zdisp = value;
	 }
	 else if (strcmp(option,"rho1") == 0) {
	    rho1 = value;
	 }
	 else if (strcmp(option,"rho2") == 0) {
	    rho2 = value;
	 }
	 else if (strcmp(option,"period_proportion") == 0) {
	    period_proportion = value;
	 }
	 else if (strcmp(option,"period_proportion_stdev") == 0) {
	    period_proportion_stdev = value;
	 }
	 else if (strcmp(option,"period_min") == 0) {
	    period_min = long(value);
	 }
	 else if (strcmp(option,"period_max") == 0) {
	    period_max = long(value);
	 }
	 else if (strcmp(option,"tilt_lattice") == 0) {
	    do_tilt_lattice = (value == 1.0);
	 }
	 else if (strcmp(option,"use_random_orientation") == 0) {
	    use_random_orientation = (value == 1.0);
	 }
	 else if (strcmp(option,"packdist_radius_fraction") == 0) {
	    packdist_radius_fraction = value;
	 }
	 else if (strcmp(option,"packdist_addition") == 0) {
	    packdist_addition = value;
	 }
	 else {
	    fprintf(stderr,"Warning: unrecognized parameter %s\n",option);
	 }

	 char temp_string[8192];
	 strcpy(temp_string,param_string);
	 sprintf(param_string,"%s%s = %f\n",temp_string,option,value);
      }
   }
   else {
      fprintf(stderr,"Warning: could not read parameters\n");
   }
}



