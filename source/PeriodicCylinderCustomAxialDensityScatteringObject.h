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

#ifndef _PERIODIC_CYLINDER_CUSTOM_AXIAL_DENSITY_SCATTERING_OBJECT_H
#define _PERIODIC_CYLINDER_CUSTOM_AXIAL_DENSITY_SCATTERING_OBJECT_H

using namespace std;

#include "ScatteringObject.h"

class PeriodicCylinderCustomAxialDensityScatteringObject : public ScatteringObject {
 public:
    ~PeriodicCylinderCustomAxialDensityScatteringObject() {};


    virtual vector<HSComplex> CalculateScatteredAmplitude(
	vector<double> &kxArray, vector<double> &kyArray,
	vector<double> &kzArray);

    virtual vector<HSComplex> CalculateScatteredAmplitudeRandomized(
	vector<double> &kxArray, vector<double> &kyArray,
	vector<double> &kzArray, long seed, long numpoints);




 private:
    PeriodicCylinderCustomAxialDensityScatteringObject() {};


    void calculate_disc_structure_factor(
	vector<HSComplex> &outStructureFactor, vector<double> &kxArray, 
	vector<double> &kyArray, vector<double> &kzArray);

    void calculate_axial_structure_factor(
	vector<HSComplex> &outStructureFactor, vector<double> &kxArray, 
	vector<double> &kyArray, vector<double> &kzArray);


    double x,y,z;
    double alpha,beta;
    double R;
    double L;
    long n;


    static void InitializeAxialScatteringAmplitude(const vector<double> &coord, 
					   const vector<double> & dens);

    // only one type of axial electron density profile can be 
    // active during the execution of the program
    static vector<HSComplex> axial_scattering_amplitude;
    static double k_step;
    const static double MAX_PERIOD_LENGTH = 1000;



 public:
    static ScatteringObject* InstantiateObject(char *s);
    

    static void SetAxialDensity(const char *filename);



};



#endif // _PERIODIC_CYLINDER_CUSTOM_AXIAL_DENSITY_SCATTERING_OBJECT_H
