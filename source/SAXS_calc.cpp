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


/*
  This program calculates scattering for simple geometrical objects
  and their arbitrary combinations. The object data and k-vector
  values are read from ascii files and the result is printed to
  another ascii file. Multiple data-files can be used, in which case
  multiple result files are also produced. 
*/

#include <unistd.h>

#include "SAXS_calc.h"
#include "ScatteringObject.h"



int main(int argc, char *argv[])
{
   strcpy(edfilename,"");
   int file_start_index = parse_options(argc,argv);

   for (int file_index=file_start_index;file_index<argc;file_index++) {
      read_data(argv[file_index]);
      vector<double> intensity(kxArray.size(),0);      

      for (long i=0;i<scatteringObjects.size();i++) {
	 vector<HSComplex> amplitude(kxArray.size());
	 
	 for (long j=0;j<scatteringObjects[i].size();j++) {
	    ScatteringObject *obj = scatteringObjects[i][j];
	    
	    vector<HSComplex> current_amplitude = obj->CalculateScatteredAmplitude(
	       kxArray,kyArray,kzArray);


	    for (long n=0;n<current_amplitude.size();n++) {
	       amplitude[n] = amplitude[n] + current_amplitude[n];
	    }
	 }

	 for (long n=0;n<amplitude.size();n++) {
	    HSComplex c = amplitude[n];
	    intensity[n] += c.Re()*c.Re() + c.Im()*c.Im();
	 }
      }


      write_result(intensity,argv[file_index]);
      clear_data();
   }

   return 0;
}


int parse_options(int argc, char *argv[])
{
   char *optstring = "k:c:";
   char option;
   bool k_read = false;

   while ((option = getopt(argc,argv,optstring)) != -1) {
      switch (option) {
      case 'k':
	{
	  char filename[1024];
	  int k = sscanf(optarg,"%s",filename);
	  
	  strcpy(kfilename,filename);

	  read_kfile(filename);
	  k_read = true;
	  break;
	}

      case 'c':
	{
	  char filename[1024];
	  int k = sscanf(optarg,"%s",filename);
	  strcpy(edfilename,filename);
	  ScatteringObject::SetCustomEDFile(filename);
	  break;
	}

      default:
	break;
      }
   }

   if (!k_read) {
      fprintf(stderr,"use -k<filename> to give file with k-values\n");
      exit(1);
   }

   return optind;
}



void read_kfile(char *filename)
{
   FILE *file = fopen(filename,"r");

   double kx,ky,kz;

   if (file != NULL) {
      while (fscanf(file,"%lf %lf %lf",&kx,&ky,&kz) == 3) {
	 kxArray.push_back(kx);
	 kyArray.push_back(ky);
	 kzArray.push_back(kz);

	 double theta = atan2(kx,sqrt(kz*kz+ky*ky));
	 double phi = atan2(kz,ky);

	 cos_theta_Array.push_back(cos(theta));
	 sin_theta_Array.push_back(sin(theta));	 

	 cos_phi_Array.push_back(cos(phi));
	 sin_phi_Array.push_back(sin(phi));
      }
   } 
   else {
      printf("can not open k-file %s\n",filename);
      exit(1);
   }

   fclose(file);
}



void read_data(char *filename)
{
   FILE *fp = fopen(filename,"r");
   char buffer[1024];
   bool inside_comment = false;


   vector<ScatteringObject*> current_objects;

   while (fgets(buffer,1024,fp) != NULL) {
      if (strncmp(buffer,"###",3) == 0) {
	 inside_comment = !inside_comment;
      }
      else if (inside_comment) {
	 strcat(comment_buffer,buffer);
      }
      else if (!inside_comment) {
	 if (buffer[0] == '*') {
	    scatteringObjects.push_back(current_objects);
	    current_objects.clear();
	 }
	 else if ((buffer[0] != EOF) && (buffer[0] != '%')) {
	    ScatteringObject *obj = ScatteringObject::Manufacture(buffer);
	
	    if (obj != NULL) {
	       current_objects.push_back(obj);
	    }
	 }
      }
   }

   if (current_objects.size() > 0) {
      scatteringObjects.push_back(current_objects);
      current_objects.clear();
   }


   fclose(fp);
}


void clear_data()
{
   strcpy(comment_buffer,"");

   for (long i=0;i<scatteringObjects.size();i++) {
      for (long j=0;j<scatteringObjects[i].size();j++) {
	 delete scatteringObjects[i][j];
	 scatteringObjects[i][j] = NULL;
      }
   }

   scatteringObjects.clear();
}


void write_result(vector<double> &intensity, char *filename)
{
   char resultfilename[1024];
   sprintf(resultfilename,"%s_SAXS_calc",filename);
   FILE *fp = fopen(resultfilename,"w+");

   if (fp != NULL) {
      write_header(fp);
      for (long i=0;i<kxArray.size();i++) {
	 fprintf(fp,"%.6f %.6f %.6f %.16f\n",
		 kxArray[i],kyArray[i],kzArray[i],intensity[i]);
      }
   }

   fclose(fp);
}


void write_header(FILE *fp)
{
   fprintf(fp,"%s",comment_buffer);
   fprintf(fp,"k file name = %s\n",kfilename);
   fprintf(fp,"electron density file name = %s\n",edfilename);
   fprintf(fp,"*********\n");
}
