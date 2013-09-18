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

#ifndef _HS_COMPLEX_H
#define _HS_COMPLEX_H

#include <math.h>
#include <stdio.h>

class HSComplex {
  public:
   HSComplex() { _re = 0; _im = 0; }
   HSComplex(const HSComplex &z) { _re = z._re; _im = z._im; }
   HSComplex(double re, double im=0) { _re = re; _im = im; }

   
   void set(double re, double im) { _re = re; _im = im; }

   HSComplex& operator=(const HSComplex &z) {
      _re = z._re; _im = z._im;
      return *this;
   }

   HSComplex& operator=(const double re) {
     _re = re; _im = 0;
   }


   HSComplex& operator+=(const HSComplex &z) {
      _re += z._re; _im += z._im;
      return *this;
   }

   HSComplex& operator-=(const HSComplex &z) {
      _re -= z._re; _im -= z._im;
      return *this;
   }

   HSComplex& operator*=(const HSComplex &z) {
      double re = _re * z._re - _im * z._im;
      double im = _re * z._im + _im * z._re;
      _re = re; _im = im;
      return *this;
   }

   HSComplex& operator-() {
      _re = -_re; _im = -_im;
      return *this;
   }

   HSComplex operator+(const HSComplex &z) const {
      HSComplex c(*this);
      c += z;
      return c;
   }
   
   HSComplex operator-(const HSComplex &z) const {
      HSComplex c(*this);
      c -= z;
      return c;
   }
   
   HSComplex operator*(const HSComplex &z) const {
      HSComplex c(*this);
      c *= z;
      return c;
   }
   

   double arg() const { return sqrt(_re*_re + _im*_im); }
   double phase() const { return atan2(_im,_re); }
   double Re() const { return _re; }
   double Im() const { return _im; }

   void print(FILE *fp) {
      fprintf(fp,"(%f,%f)\n",_re,_im);
   }

  private:
   double _re,_im;
};

#endif // _HS_COMPLEX_H






