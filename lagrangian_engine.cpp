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
#include "constants.h"
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

int RK4(double t0, vectorXYZ *point, int (*velocity)(double t,vectorXYZ point, vectorXYZ *vint));
int FlowVplusSinkV(double t,vectorXYZ point, vectorXYZ *vint);

// This function prints a little command line manual if there are invalid parameters 
void print_usage(const string me) 
{
  cout << "Usage: " << me << " <file_parameters>" << endl << endl;
}

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

  int (*velocity)(double ,vectorXYZ , vectorXYZ* );

  //velocity = GetVelocity; // choice of velocity equation = veloctiy field
  velocity = FlowVplusSinkV; // choice of velocity particle = flow velocity + sinking velocity

#pragma omp parallel for default(shared) private(t) // Parallelizing the code for computing trajectories

  for (i = 0; i < numtracers; i++)
    {
      for(t=0; t<tau-intstep; t+=intstep)
	{	
	  // Semi-implicit 4th order Runge-Kutta
	  if(RK4(t, &tracer[i], velocity)!=0|| tracer[i].z <= finaldepth)
	    {
	      timespent[i] = t;
	      break;
	    }
	  timespent[i] = t+intstep;

	  if(tracer[i].z <= finaldepth)	      
	    break;
	    
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
int RK4(double t0, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* ))
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
  if(velocity(t0,*point, &v1))
    return 1;
  h = rearth * cos(rads*(point->y));
  v1.x = degrees*(v1.x / h ); // rads velocity
  v1.y = degrees*(v1.y / rearth); // rads velocity


  /* Calculate V2: */
  t = t0 + tstep2;
  point2 = *point + (tstep2 * v1);

  if(velocity(t,point2, &v2))
    return 1;

  h = rearth * cos(rads*(point2.y));
  v2.x = degrees*(v2.x / h); // rads velocity
  v2.y = degrees*(v2.y / rearth); // rads velocity


  /* Calculate V3: */
  point3 = *point + (tstep2 * v2);

  if(velocity(t,point3, &v3))
    return 1;

  h = rearth * cos(rads*(point3.y));
  v3.x = degrees*(v3.x / h);
  v3.y = degrees*(v3.y / rearth);

  
  /* Calculate V4: */
  t = t0 + intstep;
  point4 = *point + (intstep * v3);

  if(velocity(t,point4, &v4))
    return 1;

  h = rearth * cos(rads*(point4.y));
  v4.x = degrees*(v4.x / h);
  v4.y = degrees*(v4.y / rearth);

  /* Calculate Final point */  
  *point += (tstep6 * (v1 + v4 + 2.0*(v2 + v3)));

  return 0;
}

int FlowVplusSinkV(double t,vectorXYZ point, vectorXYZ *vint)
{
  if(GetVelocity( t, point, vint))
    return 1;
  
  vint->z = vint->z - vsink;
  return 0;
}
