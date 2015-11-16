#include <iostream>
#include <iomanip>
#include <ctime>  
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

#include "parameters.h" // Function to read parameters.dat
#include "ioutil.h"
#include "constants.h"
#include "velocity.h" // Function to read velocities 
#include "VTK.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
// Declaring constants and extern variables
/////////////////////////////////////////////////////////////////////////////////////////////////



extern int verbose;
extern double nodesize;
extern double particlespacing;
extern double startdepth;
extern double finaldepth;
extern date reference_date;

extern double network_ll_lat;
extern double network_ll_lon;
extern double network_tr_lat;
extern double network_tr_lon;

// these parameters are used to name files. Better do this once time in paramters.cpp and use the name files directly 
extern char velocitydir[];
extern string filename_upnetwork;
extern string filename_downnetwork;
extern string filename_itracer;

class  rectangle 
{
 public:
  vectorXYZ center;
  vectorXYZ tr;
  vectorXYZ ll;
  double land_ratio;
};

vector<rectangle> upnode;
vector<rectangle> downnode;

// total number of nodes
int num_nodes=0; 

// Variables
// coordinates of points of the grid
double points_long_vec;
double points_lat_vec;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN CODE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function prints a little command line manual if there are invalid parameters 
void print_usage(const string me) 
{
  cout << "Usage: " << me << " <file_parameters>" << endl << endl;
}


int main( int argc, char *argv[] )
{
  double lon,lat;
  double latmin,latmax;
  double lonmin,lonmax;
  double delta_lon, delta_lat;

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
  // Read domain parameters from parameters.dat (grid step, low-left long, low-left lat, top-right lon, top-right lat)

  if(readparams(fnameparameters)==1)
    {
      // If it is not specified the path to parameters file
      cout << "Error in reading file parameters" << endl;// Show a error message 
      return 1;
    }

 // Reading the velocity field and longitude and latitude in the velocity file correponds to startdate
  if(LoadVelocityGrid(reference_date,velocitydir)!=0)
    cout << "Error in reading velocities"<< endl;
  
  /* VTK OUTPUT */
  string name_VTKfile_bathymetry="bathymetry.vtk";

  MakeVTKStructuredGrid2D(vgrid_lon, vgrid_lat, startdepth, bathymetry, nlon, nlat, name_VTKfile_bathymetry );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // CONSTRUCTION of NETWORK GRID                                                                               //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(verbose == 1) // Print screen for debugging
    {
      cout <<"GRID CONSTRUCTION PARAMETERS:"<<endl; 
      cout << " Grid step is:"<< nodesize<<endl;
      cout << " Low left grid corner (lat,long)=(" << network_ll_lat << "," << network_ll_lon << ")"<<endl;
      cout << " Top right grid corner (lat,long)=(" << network_tr_lat << "," << network_tr_lon << ")"<<endl;
    }


  ofstream ofile_upnetwork;
  // opening output files
  if(FileExist(filename_upnetwork.c_str()))
    {
      cout << "el archivo existe"<<endl;
    }  
  ofile_upnetwork.open(filename_upnetwork.c_str());

  ////////////////////////////////////////////////////////////////////////
  ////////// Computing center of nodes with 2 loops ////////////////////// 
  ////////////////////////////////////////////////////////////////////////
  int i_ll;
  int j_ll;
  int i_tr;
  int j_tr;
  int land = 0;
  int point = 0;

  num_nodes = 0;

  int capacity_node;
  
  capacity_node = ((network_tr_lat-network_ll_lat)*(network_tr_lon-network_ll_lon)) /(nodesize*nodesize);

  upnode.reserve(capacity_node);//reserve space for allocation, it reduces the execution time

  // Top right corner of network grid = (network_tr_lat,network_tr_lon)
  // Low left corner of network grid = (network_ll_lat,network_ll_lon)
  // Top right center of network grid = (latmax,lonmax)
  // Low left center of network grid = (latmin,lonmin) 
  delta_lat = nodesize;
  latmin = network_ll_lat + (delta_lat/2.0);
  latmax = network_tr_lat - (delta_lat/2.0);
  for(lat=latmin; lat<latmax; lat+=delta_lat)
    {
      delta_lon=nodesize/(cos(rads*lat));
      lonmin = network_ll_lon + (delta_lon/2.0);
      lonmax = network_tr_lon - (delta_lon/2.0);
      
      for(lon=lonmin; lon<lonmax; lon+=delta_lon)
	{
	  
	  // Initialization
	  land=0;
	  point=0;

	  ////////// Computing ratio of land/sea for each node /////////////
	  upnode.push_back(rectangle());
	  upnode[num_nodes].center.x = lon;
	  upnode[num_nodes].center.y = lat;
	  upnode[num_nodes].center.z = startdepth;

	  // finding two corners of the node from the center
	  upnode[num_nodes].ll.x = lon - delta_lon/2.0;
	  upnode[num_nodes].ll.y = lat - delta_lat/2.0;
	  upnode[num_nodes].ll.z = startdepth;

	  upnode[num_nodes].tr.x = lon + delta_lon/2.0;
	  upnode[num_nodes].tr.y = lat + delta_lat/2.0;
	  upnode[num_nodes].tr.z = startdepth;
	  
	  // co-locate the coordinate of the node in the grid of the velocity field  
	  // Lower left corner, outside
	  i_ll = GetLonIndex(upnode[num_nodes].ll.x);
	  j_ll = GetLatIndex(upnode[num_nodes].ll.y);
	  // Top right, inside
	  i_tr = GetLonIndex(upnode[num_nodes].tr.x);
	  j_tr = GetLatIndex(upnode[num_nodes].tr.y);
	  // Top right, outside
	  i_tr = i_tr + 1;
	  j_tr = j_tr + 1;
	  
	  // count in velocity field the number of missing & real values 
	  for(int i=i_ll; i<=i_tr; i+=1)
	    {
	      for(int j=j_ll; j<=j_tr; j+=1)
		{

		  //test if missing value (1.e+20 = missing values for MyOcean)
		  if(startdepth < bathymetry[i][j])
		    { 
		      land++;
		    }
		  point++;
		}
	    }
	  
	  // compute proportion of land per node (land_ratio=1 means land only)
	  upnode[num_nodes].land_ratio = ((double) land)/((double) point);
			    
	  // For each node (1 node = 1 line in output file): write lon, lat of node center + longitude extension + land ratio
	  ofile_upnetwork << upnode[num_nodes].center.x << " " << upnode[num_nodes].center.y << " " << delta_lon << " " << upnode[num_nodes].land_ratio << endl;
	  num_nodes++;
	}
    }

  ofile_upnetwork.close();

  // Print out total number of nodes if verbose = 1 
  if(verbose == 1)
    {
      cout <<" Total number of nodes: "<< num_nodes << endl;   
      cout <<" Total number of nodes <vector>: "<< upnode.size() << endl;   
      cout <<" Capacity nodes <vector>: "<< upnode.capacity() << endl;   
    }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  CONSTRUCTION of FINAL NETWORK          
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // opening output files
  ofstream ofile_downnetwork;
  ofile_downnetwork.open(filename_downnetwork.c_str());

  int q;
  downnode.reserve(capacity_node);//reserve space for allocation, it reduces the execution time
  for(q=0; q<num_nodes; q++)
    {
	  
      // Initialization
      land=0;
      point=0;
      
      ////////// Computing ratio of land/sea for each node /////////////
      downnode.push_back(rectangle());
      downnode[q] = upnode[q];
      
      upnode[num_nodes].center.z = finaldepth;
      upnode[num_nodes].ll.z = finaldepth;
      upnode[num_nodes].tr.z = finaldepth;
      
      // co-locate the coordinate of the node in the grid of the velocity field  
      // Lower left corner, outside
      i_ll = GetLonIndex(downnode[q].ll.x);
      j_ll = GetLatIndex(downnode[q].ll.y);
      // Top right, inside
      i_tr = GetLonIndex(downnode[q].tr.x);
      j_tr = GetLatIndex(downnode[q].tr.y);
      // Top right, outside
      i_tr = i_tr + 1;
      j_tr = j_tr + 1;
      
	  // count in velocity field the number of missing & real values 
	  for(int i=i_ll; i<=i_tr; i+=1)
	    {
	      for(int j=j_ll; j<=j_tr; j+=1)
		{

		  //test if missing value (1.e+20 = missing values for MyOcean)
		  if(finaldepth < bathymetry[i][j])
		    { 
		      land++;
		    }
		  point++;
		}
	    }
	  
	  // compute proportion of land per node (land_ratio=1 means land only)
	  downnode[q].land_ratio = ((double) land)/((double) point);
			    
	  // For each node (1 node = 1 line in output file): write lon, lat of node center + longitude extension + land ratio
	  ofile_downnetwork << downnode[q].center.x << " " << downnode[q].center.y << " " <<  downnode[q].tr.x-downnode[q].ll.x << " " << downnode[q].land_ratio << endl;
    }

  ofile_downnetwork.close();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////  CONSTRUCTION of INITIAL POSITIONS OF TRACERS //////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // opening output files
  ofstream ofile_itracer;
  ofile_itracer.open(filename_itracer.c_str());

  ////////////////////////////////////////////////////////////////////////
  ////// Computing tracer initial position with 2 loops //////////////////
  ////////////////////////////////////////////////////////////////////////

  // Start/end seeding tracers
  latmin = network_ll_lat;
  latmax = network_tr_lat;
  delta_lat = particlespacing;

  lonmin = network_ll_lon;
  lonmax = network_tr_lon;

  if(verbose == 1) // Print screen for debugging
    {
      cout <<"TRACERS CONSTRUCTION PARAMETERS:"<<endl; 
      cout << " Particle latspacing:"<< particlespacing <<endl;
      cout << " Low left grid corner (lat,long)=(" << latmin << "," << lonmin << ")"<<endl;
      cout << " Top right grid corner (lat,long)=(" << latmax << "," << lonmax << ")"<<endl;
    }

  int num_tracers = 0;
  for(lat=latmin; lat<latmax; lat+=delta_lat)
    {
      
      delta_lon=particlespacing/(cos(rads*lat));         

      for(lon=lonmin; lon<lonmax; lon+=delta_lon)
	{

	  // Find the 4 nearest neighboors of the tracer in the grid of the velocity field  
	  // Lower left corner 
	  i_ll = GetLonIndex(lon);
	  j_ll = GetLatIndex(lat);

	  if((i_ll >= nlon-1) || (j_ll >=nlat-1))
	    continue;

	  // Top right corner
	  i_tr = i_ll + 1;
	  j_tr = j_ll + 1;
	  
	  // initialize tracer only if there at least 1 or more nearest neighboors with velocity (no initialization when 4 neighboors = NaN) 
	  if((startdepth >= bathymetry[i_ll][j_ll]) || 
	     (startdepth >= bathymetry[i_ll][j_tr]) || 
	     (startdepth >= bathymetry[i_tr][j_ll]) || 
	     (startdepth >= bathymetry[i_tr][j_tr]))
	    {
	      // Print out initial positions (lat, lon) of tracers (1 line per particle in output files) 
	      ofile_itracer << lon << " " << lat << endl;
	      num_tracers++;
	    }
	}
    }
  ofile_itracer.close();

  // Print out total number of particles if verbose = 1 
  if(verbose==1)
    {
      cout <<" Total number of tracers: "<< num_tracers << endl;
    } 

  FreeMemoryVelocityGrid();

  return 0;
}

