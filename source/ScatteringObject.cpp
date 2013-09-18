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

#include "PeriodicCylinderScatteringObject.h"
#include "PeriodicCylinderCustomAxialDensityScatteringObject.h"
#include "ScatteringObject.h"


ScatteringObject* ScatteringObject::Manufacture(char *s)
{
   char type[32];

   int n;
   sscanf(s,"%s %n",type,&n);

   ScatteringObject *o = NULL;
   
   if (strcasecmp("cylper",type) == 0) {
      o = PeriodicCylinderScatteringObject::InstantiateObject(s+n);
   }
   else if(strcasecmp("cylpercad",type) == 0) {
     o = PeriodicCylinderCustomAxialDensityScatteringObject::InstantiateObject(s+n);
   }
   else {
      fprintf(stderr,"Warning: Unknown object type %s\n",type);
   }

   return o;
}



void ScatteringObject::SetCustomEDFile(const char *filename)
{
  PeriodicCylinderCustomAxialDensityScatteringObject::SetAxialDensity(filename);
}
