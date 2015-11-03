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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////// prototype-functions used in the main code (written below, at the end of the main code)
/////////////////////////////////////////////////////////////////////////////////////////////////

string DoubletoInt(int ndigits, int ndecimals, double number);
int SinkRK4(double t0, vectorXYZ *point);

// This function prints a little command line manual if there are invalid parameters 
void print_usage(const string me) 
{
  cout << "Usage: " << me << " <file_parameters>" << endl << endl;
}

// This function counts the number of lines of a file
int countlines(string namefile)
  {
    ifstream ifile(namefile.c_str());
    int numlines;
    string line;
    
    if (!ifile.is_open()) 
      {
	string me = "countlines()";
	cout << me <<": Skipping unreadable file \"" << namefile.c_str() << "\" "<<endl;
	return 1;
      }
    for ( numlines = 0; getline(ifile, line); ++numlines)
      ;

    ifile.close();
    return numlines;
  }

// Directories Outputs
//char const output_dir[] = "./";
//char const input_dir[] = "./";

//Macros
#define TRIAL_POINT(rtrial, r, alpha, V)	  \
  rtrial.x = r.x + (alpha) * V.x;		  \
  rtrial.y = r.y + (alpha) * V.y;		  \
  rtrial.z = r.z + (alpha) * V.z

#define R_EARTH 6371000.0

/////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////
// Declaring constant and global variables
/////////////////////////////////////////////////////////////////////////////////////////////////

// Conversions
//const double day2sec = 60. * 60. * 24.; // Daily data -> velocities converted from seconds to days
//const double m2cm = 100.; // Input velocities in cm -> m for the interpolation
//const double pig = 3.14159265358;  // pi greca
//const double deg2rad = (pig/180.); // Conversion from degree to radians
//const double Radius = 6372795.48;     // radius of the Earth (m)

/////////////////////////////////////////////////////////////////////////////////////////////////

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

  string sdomain_network;
  string svelocity_field;
  string sdepth;
  string sparticle_latspacing;

  string sstart_date;
  string stau;
  string sint_step;
  string name_input_tracer;
  string name_output_tracer;
  char buffer[50];

  // reconstruct name of the input file (tracer initial positions)
  sprintf(buffer, "networkdomain%d", domain_network);
  sdomain_network = buffer;

  sprintf(buffer, "_vflow%d", velocity_field);
  svelocity_field=buffer;

  sdepth="_depth";
  sparticle_latspacing = "_particlelatspacing";

  name_input_tracer = 
    input_dir + 
    sdomain_network +
    svelocity_field + 
    sdepth + DoubletoInt(4, 0, start_depth) + 
    sparticle_latspacing + DoubletoInt(5,4,particle_latspacing) + 
    postfixitracer; 

  if(verbose == 1) // Print input file name (tracer initial positions) to screen for debugging
    {  
      cout << " Input file name (tracer initial positions) "<< name_input_tracer << endl;
    }

  // Call function to read the number of lines in input file name .trac (tracer initial positions)
  int numtracers=countlines(name_input_tracer);
  if(numtracers < 0)
    {
      cout << "ERROR counting lines of file "<< name_input_tracer << endl;
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
  ifstream file_tracers(name_input_tracer.c_str());

  int i;
  
  // Read line by line the lon/lat for all initial positions of all particles 
  for(i=0; i<numtracers; i++)
    {
      file_tracers >> tracer[i].x >> tracer[i].y;
      tracer[i].z = start_depth;
    }
  file_tracers.close();

  // READ layer numbered "layer_index" of the VELOCITY FIELD (from initial day "start_date" to "start_date + tau + 1")

  if(LoadVelocityGrid(reference_date)!=0)
    {
      cout << "Error in reading reference date netcdf"<< endl;
      return 1;
    }

  if(LoadVelocities(start_date, tau)!=0)
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
      for(t=0; t<tau-1; t+=int_step)
	{	
	  // Semi-implicit 4th order Runge-Kutta
	  if(SinkRK4(t, &tracer[i])!=0|| tracer[i].z <= final_depth)
	    {
	      timespent[i] = t-int_step;
	      break;
	    }
	  timespent[i] = t;
	}
    }

  // Constructing name of output file (= final positions of tracers)
  
  sstart_date="_startdate";

  sint_step = "_intstep";

  string sidepth = "_idepth";
  string sfdepth = "_fdepth";
  string svsinking = "_vsinking";

  name_output_tracer = 
    output_dir + 
    sdomain_network +
    svelocity_field + 
    sidepth + DoubletoInt(4, 0, start_depth) + 
    sfdepth + DoubletoInt(4, 0, final_depth) + 
    sparticle_latspacing + DoubletoInt(5,4,particle_latspacing) +
    sstart_date + DoubletoInt(2,0,start_date.day) + DoubletoInt(2,0,start_date.month) + DoubletoInt(4,0,start_date.year) +
    svsinking + DoubletoInt(3,0,v_sinking) +
    sint_step +  DoubletoInt(4,2,int_step) +
    postfixftracer;
  
  // open output file (tracer final positions) 
  ofstream file_outtracers(name_output_tracer.c_str());// open tracer file 

  // write final positions
  for (i = 0; i < numtracers; i++)
    {
      file_outtracers << tracer[i].x <<" "<< tracer[i].y<<" "<< tracer[i].z<<" "<< timespent[i] <<endl; 
    }
  file_outtracers.close();
  
  // Free memory
  delete [] tracer;
  delete [] timespent;
  FreeMemoryVelocityGrid();
  FreeMemoryVelocities();
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
  tstep2 = int_step*0.5;
  tstep6 = int_step/6.0;
  
  /* Calculate V1: */
  if(GetVelocity(t0,*point, &v1))
    return 1;
  h = R_EARTH * cos(RADS((*point).y));
  v1.x = DEGREE(v1.x / h ); // rads velocity
  v1.y = DEGREE(v1.y / R_EARTH); // rads velocity
  v1.z = v1.z - v_sinking; //adding sinking velocity

  /* Calculate V2: */
  t = t0 + tstep2;
  TRIAL_POINT(point2, (*point), tstep2, v1);
  if(GetVelocity(t,point2, &v2))
    return 1;

  h = R_EARTH * cos(RADS(point2.y));
  v2.x = DEGREE(v2.x / h); // rads velocity
  v2.y = DEGREE(v2.y / R_EARTH); // rads velocity
  v2.z = v2.z - v_sinking;//adding sinking velocity

  /* Calculate V3: */
  TRIAL_POINT(point3, (*point), tstep2, v2);
  if(GetVelocity(t,point3, &v3))
    return 1;

  h = R_EARTH * cos(RADS(point3.y));
  v3.x = DEGREE(v3.x / h);
  v3.y = DEGREE(v3.y / R_EARTH);
  v3.z = v3.z - v_sinking;//adding sinking velocity
  
  /* Calculate V4: */
  t = t0 + int_step;
  TRIAL_POINT(point4, (*point), int_step, v3);
  if(GetVelocity(t,point4, &v4))
    return 1;

  h = R_EARTH * cos(RADS(point4.y));
  v4.x = DEGREE(v4.x / h);
  v4.y = DEGREE(v4.y / R_EARTH);
  v4.z = v4.z - v_sinking;

  /* Calculate Final point */
  point->x = point->x + tstep6 * (v1.x + v4.x + 2.0 * (v2.x + v3.x));
  point->y = point->y + tstep6 * (v1.y + v4.y + 2.0 * (v2.y + v3.y));
  point->z = point->z + tstep6 * (v1.z + v4.z + 2.0 * (v2.x + v3.z));

  return 0;
}
string DoubletoInt(int ndigits, int ndecimals, double number)
{
  ostringstream stream;// it needs to include  <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>

  stream << fixed; //it needs to include <iostream>

  stream << setfill('0') << setw(ndigits);

  stream << setprecision(0) << (abs(number)*factor);
  
  return stream.str();
}

