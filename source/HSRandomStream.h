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

#ifndef _HS_RANDOM_STREAM_H
#define _HS_RANDOM_STREAM_H

#include <math.h>


class HSRandomStream {
public:
    HSRandomStream(long seed) {
	_prev_normal_available = false;
	Initialize(seed);
    }

    inline void Initialize(long seed);



    inline long discrete_uniform(long a,long b); // returns a uniformly selected integer in [a,b] (must be a <= b)
    inline double continuous_uniform(double a, double b); // returns a uniformly selected value in ]a,b[


    inline double normal_distribution(); // returns a value from distribution 1/sqrt(pi)*exp(-x^2)
    inline double gaussian_distribution(double mu, double stdev); // returns a value from distribution exp(-(x-mu)^2/stdev^2)
    inline double exponential_distribution(double lambda); // returns a value from distribution lambda*exp(-lambda*x)

    inline double lognormal_distribution(double mean, double var);
    inline double triangular_distribution(double a, double b, double c); // a <= b <= c
    inline double weibull_distribution(double alpha, double beta); // alpha > 0 (shape), beta > 0 (scale)

    inline double rand1(); // returns uniform value in ]0,1[



private:
    static const long NTAB = 32;
    long _iy1;
    long _iv_table[NTAB];
    long _idum;

    static const long IA1 = 16807, IM1 = 2147483647, IQ1 = 127773, IR1 = 2836;
    static const long NDIV1 = (1 + (IM1 - 1)/NTAB);
    static const double EPS = 3.0e-16,AM1 = 1.0 / IM1, RNMX = (1.0 - EPS);


    bool _prev_normal_available;
    double _prev_normal;
   

};




// this returns an integer in [a,b] that is uniformly selected
long HSRandomStream::discrete_uniform(long a, long b)
{
   double x = rand1();

   return long(a + floor((b-a+1) * x));
}



// this returns a uniformly selected value in ]a,b[
double HSRandomStream::continuous_uniform(double a, double b)
{
   return a + (b-a) * rand1();
}




double HSRandomStream::normal_distribution()
{
   if (_prev_normal_available) {
      _prev_normal_available = false;
      return _prev_normal;
   }
   else {
      double x1,x2;
      double rsq;
      do {
	 x1 = 2.0 * rand1() - 1.0;
	 x2 = 2.0 * rand1() - 1.0;

	 rsq = x1 * x1 +  x2 * x2;
      } while ((rsq >= 1) || (rsq == 0));

      double fac = sqrt(-2*log(rsq)/rsq);

      _prev_normal = fac * x1;
      _prev_normal_available = true;

      return fac * x2;
   }
}


double HSRandomStream::gaussian_distribution(double mu, double a)
{
   double x = normal_distribution();

   return mu + x * a;
}


double HSRandomStream::exponential_distribution(double lambda)
{
   double x = rand1();

   while (x == 0)
      x = rand1();

   return -log(x) / lambda;
}



double HSRandomStream::lognormal_distribution(double mean, double var)
{
   double m1 = mean;
   double v1 = var;
   double m = log(m1*m1 / sqrt(m1*m1+v1));
   double v = log((m1*m1 + v1)/(m1*m1));

   double Y = gaussian_distribution(m,sqrt(v));
   return exp(Y);
}



double HSRandomStream::triangular_distribution(double a, double b, double c)
{
   double x = rand1();

   double bp = (b-a) / (c-a);


   if (x < bp) {
      x = sqrt(x * bp);
   }
   else {
      x = 1 - sqrt((1-bp)*(1-x));
   }

   return a+(c-a)*x;
}


double HSRandomStream::weibull_distribution(double alpha, double beta)
{
   double u = rand1();

   return beta * pow(-log(u),1/alpha);
}


// this is basically ran1 from Numerical Recipes in C++ p. 284
// returns a value in the range ]0,1.0[
double HSRandomStream::rand1()
{
   double temp;
   long j,k = _idum / IQ1;
   _idum = IA1 * (_idum - k * IQ1) - IR1 * k;
   if (_idum < 0) 
      _idum += IM1;

   j = _iy1 / NDIV1;
   _iy1 = _iv_table[j];
   _iv_table[j] = _idum;

   temp = AM1 * _iy1;
   if (temp > RNMX)
      return RNMX;
   else
      return temp;
}




void HSRandomStream::Initialize(long seed) 
{
   // initialize for rand1
   _idum = (seed > 0 ? seed : 1);
   
   for (long j=NTAB+7;j>=0;j--) {
      long k = _idum / IQ1;
      _idum = IA1 * (_idum - k * IQ1) - IR1 * k;
      if (_idum < 0) 
	 _idum += IM1;
      if (j < NTAB)
	 _iv_table[j] = _idum;
   }
   _iy1 = _iv_table[0];
}






#endif // _HS_RANDOM_STREAM_H
