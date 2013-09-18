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
#ifndef _HS_3D_VECTOR_H
#define _HS_3D_VECTOR_H

#include <math.h>
#include <stdio.h>

class HS3DVector {
   public:
      HS3DVector(double x, double y, double z) {
	 _x = x; _y = y; _z = z;
      };

      HS3DVector(const HS3DVector &v) {
	 _x = v._x; _y = v._y; _z = v._z;
      }

      HS3DVector& operator=(const HS3DVector &v) {
	 _x = v._x; _y = v._y; _z = v._z;
	 return *this;
      }

      HS3DVector& operator+=(const HS3DVector &v) {
	 _x += v._x; _y += v._y; _z += v._z;
	 return *this;
      }

      HS3DVector& operator-=(const HS3DVector &v) {
	 _x -= v._x; _y -= v._y; _z -= v._z;
	 return *this;
      }

      HS3DVector& operator*=(double a) {
	 _x = a*_x; _y = a*_y; _z = a*_z;
	 return *this;
      }


      HS3DVector operator+(const HS3DVector &v) const {
	 HS3DVector b(*this);
	 b += v;
	 return b;
      }

      HS3DVector operator-(const HS3DVector &v) const {
	 HS3DVector b(*this);
	 b -= v;
	 return b;
      }

      HS3DVector operator-() const {
	 return HS3DVector(-_x,-_y,-_z);
      }

      HS3DVector operator*(double a) const {
	 HS3DVector b(*this);
	 b *= a;
	 return b;
      }


      bool operator<(const HS3DVector &v) const {
	return (length() < v.length());
      }
     

      HS3DVector cross_product(const HS3DVector &v) const {
	 return HS3DVector(_y*v._z - _z*v._y, 
			   _z*v._x - _x*v._z, 
			   _x*v._y - _y*v._x);
      }

      double dot_product(const HS3DVector &v) const {
	 return _x*v._x + _y*v._y + _z*v._z;
      }
   
      double length() const {
	 return sqrt(_x*_x + _y*_y + _z*_z);
      }

      double dist(const HS3DVector &v) const {
	double x = _x-v._x;
	double y = _y-v._y;
	double z = _z-v._z;

	return sqrt(x*x+y*y+z*z);
      }

      double X() const {
	 return _x;
      }
      
      double Y() const {
	 return _y;
      }

      double Z() const {
	 return _z;
      }

      void print(FILE *fp) const {
	 fprintf(fp,"%.10f %.10f %.10f\n",_x,_y,_z);
      }

   private:
      double _x,_y,_z;
};



#endif // _HS_3D_VECTOR_H
