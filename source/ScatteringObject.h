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


#ifndef _SCATTERING_OBJECT_H
#define _SCATTERING_OBJECT_H


/*
  This is a parent class for scattering objects. 

  The subclasses must implement a static function that constructs the object 
  from a string and returns a pointer to this object (or NULL if the string 
  is invalid). The name of this function is InstantiateObject.

  The subclass must implement the purely virtual function 
  CalculateScatteredAmplitude that calculates the scattered amplitude for
  given values of k-vector.

  The parent class implements a static function name Manufacture 
  that takes a string that specifies the scattering object and returns 
  a pointer to the  object that was instantiated.
*/


#include <vector>
#include "HSComplex.h"


using namespace std;

class ScatteringObject {
 public:
    ScatteringObject() {};
    virtual ~ScatteringObject() {};

    static ScatteringObject* Manufacture(char *s);
    static void SetCustomEDFile(const char *filename);

    virtual vector<HSComplex> CalculateScatteredAmplitude(
	vector<double> &kxArray, vector<double> &kyArray,
	vector<double> &kzArray) = 0;

      
    virtual vector<HSComplex> CalculateScatteredAmplitudeRandomized(
	vector<double> &kxArray, vector<double> &kyArray,
	vector<double> &kzArray, long seed, long numpoints) {}; // Not implemented in the superclass


    static ScatteringObject* InstantiateObject(char *s) {return NULL;} 
};


#endif // _SCATTERING_OBJECT_H
