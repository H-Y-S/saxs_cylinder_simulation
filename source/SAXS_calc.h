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

#ifndef _SAXS_CALC_H
#define _SAXS_CALC_H

#include <vector>

#include "ScatteringObject.h"



char kfilename[1024];
char edfilename[1024];
char comment_buffer[65536];


int parse_options(int argc, char *argv[]);
void read_kfile(char *filename);
void read_data(char *filename);
void clear_data();
void write_result(vector<double> &intensity, char *filename);
void write_header(FILE *fp);

vector<double> kxArray;
vector<double> kyArray;
vector<double> kzArray;

vector<double> sin_theta_Array;
vector<double> cos_theta_Array;
vector<double> sin_phi_Array;
vector<double> cos_phi_Array;


vector<vector<ScatteringObject*> > scatteringObjects;


#endif // _SAXS_CALC_H
