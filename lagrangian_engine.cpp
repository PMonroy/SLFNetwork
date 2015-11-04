#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>

using namespace std;
#include "parameters.h" // Function to read parameters.dat
#include "velocity.h" // Function to read velocities 
#include "ioutil.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Declaring constants and extern variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern int verbose;

extern double startdepth;
extern double finaldepth;
extern date startdate;
extern double intstep;
extern int tau;
extern date reference_date;
extern double vsink;

extern char velocitydir[];
extern string filename_itracer;
extern string filename_ftracer;


/////////////////////////////////////////////////////////////////////////////////////////////////
/////// prototype-functions used in the main code (written below, at the end of the main code)
/////////////////////////////////////////////////////////////////////////////////////////////////

int SinkRK4(double t0, vectorXYZ *point);

// This function prints a little command line manual if there are invalid parameters 
void print_usage(const string me) 
{
  cout << "Usage: " << me << " <file_parameters>" << endl << endl;
}

//Macros
#define TRIAL_POINT(rtrial, r, alpha, V)	  \
  rtrial.x = r.x + (alpha) * V.x;		  \
  rtrial.y = r.y + (alpha) * V.y;		  \
  rtrial.z = r.z + (alpha) * V.z

#define R_EARTH 6371000.0
/////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN CODE 
/////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
  // READ COMMAND LINE ARGUMENTS:
  // Usage in command line when executing: ./(path executable) (path to parameters file) 
  // (ex. ./grid_construction.out parameters.dat)
  
  if (2 != argc) 
    {
      // If it is not specified the path to parameters file
      cout << "Parameter file not specified" << endl; // Show a error message  
      string me=argv[0];
      print_usage(me); // and a usage reminder (me = name of executable) 
      return 1;
    }

  // READ PARAMETERS FROM FILE
  string fnameparameters = argv[1]; // Path to parameters file = Second argument in the command line 
  if(readparams(fnameparameters)==1)
    {
      // If it is not specified the path to parameters file
      cout << "Error in reading file parameters "<< fnameparameters << endl;// Show a error message
      return 1;
    }
  
  //READ INITIAL POSITIONS

  if(verbose == 1) // Print input file name (tracer initial positions) to screen for debugging
    {  
      cout << " Input file name (tracer initial positions) "<< filename_itracer << endl;
    }

  // Call function to read the number of lines in input file name .trac (tracer initial positions)
  int numtracers=CountLines(filename_itracer.c_str());
  if(numtracers < 0)
    {
      cout << "ERROR counting lines of file "<< filename_itracer << endl;
      return 1;
    }

  if(verbose == 1) // Print total number of tracers to screen for debugging
    {
      cout <<" Total Num of particles = "<< numtracers << endl;
    }
  
  vectorXYZ *tracer;
  double *timespent;

  // Initialize vectors of lon and lat for all initial positions of all particles 
  tracer = new vectorXYZ [numtracers];
  timespent = new double [numtracers];

  // Open input file .trac (tracer initial positions)
  ifstream ifile_itracer(filename_itracer.c_str());

  int i;
  
  // Read line by line the lon/lat for all initial positions of all particles 
  for(i=0; i<numtracers; i++)
    {
      ifile_itracer >> tracer[i].x >> tracer[i].y;
      tracer[i].z = startdepth;
    }
  ifile_itracer.close();

  // READ layer numbered "layer_index" of the VELOCITY FIELD (from initial day "startdate" to "startdate + tau + 1")

  if(LoadVelocityGrid(reference_date, velocitydir)!=0)
    {
      cout << "Error in reading reference date netcdf"<< endl;
      return 1;
    }

  if(LoadVelocities(startdate, tau, velocitydir)!=0)
    {
      cout << "Error in reading velocities"<< endl;
      return 1;
    }
  
  //////////////////////////////////////////////////////////////////
  /////////////////// MAIN EVOLUTION LOOP //////////////////////////
  //////////////////////////////////////////////////////////////////

  if(verbose == 1) // Print msg to screen for debugging
    {
      cout << "Starting integration with opnemp...\n";
    }
  
  
  double t;

  //////////////////////////////////////

  // Parallelizing the code for computing trajectories 
#pragma omp parallel for default(shared) private(t)
for (i = 0; i < numtracers; i++)
    {
      for(t=0; t<tau-1; t+=intstep)
	{	
	  // Semi-implicit 4th order Runge-Kutta
	  if(SinkRK4(t, &tracer[i])!=0|| tracer[i].z <= finaldepth)
	    {
	      timespent[i] = t-intstep;
	      break;
	    }
	  timespent[i] = t;
	}
    }

  // open output file (tracer final positions) 
  ofstream ofile_ftracer(filename_ftracer.c_str());// open tracer file 
  // write final positions
  for (i = 0; i < numtracers; i++)
    {
      ofile_ftracer << tracer[i].x <<" "<< tracer[i].y<<" "<< tracer[i].z<<" "<< timespent[i] <<endl; 
    }
  ofile_ftracer.close();
  
  // Free memory
  delete [] tracer;
  delete [] timespent;
  FreeMemoryVelocityGrid();
  FreeMemoryVelocities(tau);
  return 0;
}
int SinkRK4(double t0, vectorXYZ *point)
{
  vectorXYZ point2,point3,point4;
  vectorXYZ v1,v2,v3,v4;
  
  double tstep2,tstep6;
  double h; // scale factor of equally spaced cordinate system
  double t;

  /* Time increments */
  tstep2 = intstep*0.5;
  tstep6 = intstep/6.0;
  
  /* Calculate V1: */
  if(GetVelocity(t0,*point, &v1))
    return 1;
  h = R_EARTH * cos(RADS((*point).y));
  v1.x = DEGREE(v1.x / h ); // rads velocity
  v1.y = DEGREE(v1.y / R_EARTH); // rads velocity
  v1.z = v1.z - vsink; //adding sinking velocity

  /* Calculate V2: */
  t = t0 + tstep2;
  TRIAL_POINT(point2, (*point), tstep2, v1);
  if(GetVelocity(t,point2, &v2))
    return 1;

  h = R_EARTH * cos(RADS(point2.y));
  v2.x = DEGREE(v2.x / h); // rads velocity
  v2.y = DEGREE(v2.y / R_EARTH); // rads velocity
  v2.z = v2.z - vsink;//adding sinking velocity

  /* Calculate V3: */
  TRIAL_POINT(point3, (*point), tstep2, v2);
  if(GetVelocity(t,point3, &v3))
    return 1;

  h = R_EARTH * cos(RADS(point3.y));
  v3.x = DEGREE(v3.x / h);
  v3.y = DEGREE(v3.y / R_EARTH);
  v3.z = v3.z - vsink;//adding sinking velocity
  
  /* Calculate V4: */
  t = t0 + intstep;
  TRIAL_POINT(point4, (*point), intstep, v3);
  if(GetVelocity(t,point4, &v4))
    return 1;

  h = R_EARTH * cos(RADS(point4.y));
  v4.x = DEGREE(v4.x / h);
  v4.y = DEGREE(v4.y / R_EARTH);
  v4.z = v4.z - vsink;

  /* Calculate Final point */
  point->x = point->x + tstep6 * (v1.x + v4.x + 2.0 * (v2.x + v3.x));
  point->y = point->y + tstep6 * (v1.y + v4.y + 2.0 * (v2.y + v3.y));
  point->z = point->z + tstep6 * (v1.z + v4.z + 2.0 * (v2.x + v3.z));

  return 0;
}
