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


#ifndef _UTILITIES_H
#define _UTILITIES_H

#include <math.h>
#include <string.h>

#ifndef max_c
#define	max_c(a,b)	(((a) > (b)) ? (a) : (b))
#endif	// max_c

#ifndef min_c
#define	min_c(a,b)	(((a) < (b)) ? (a) : (b))
#endif	// min_c


#ifndef round
#define round(a)	(floor((a) + 0.5))
#endif	// round



struct ltstr {  
	bool operator()(const char* s1, const char* s2) const  {
    	return strcmp(s1, s2) < 0;  
    }
};




#endif // _UTILITIES_H

