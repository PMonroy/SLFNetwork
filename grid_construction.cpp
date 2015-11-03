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
#include "velocity.h" // Function to read velocities 
#include "VTK.h"

// structures

struct rectangle {
  struct vectorXYZ center;
  struct vectorXYZ tr;
  struct vectorXYZ ll;
  double land_ratio;
};

vector<rectangle> upnode;
vector<rectangle> downnode;

// Parameters
//const double pig = 3.14159265358;
//const double pig180 = (pig/180.);

// Directories Outputs
//char const output_dir[] = "./";

// File names Outputs
//const string postfixgrid = ".grid";
//const string postfixtracer = ".trac";

// total number of nodes
int num_nodes=0; 

// Variables
// coordinates of points of the grid
double points_long_vec;
double points_lat_vec;

//Functions
string DoubletoInt(int ndigits, int ndecimals, double number);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN CODE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  string sdomain_network;
  string snum_nodes;
  string snum_tracers;
  string svelocity_field;
  string sdepth;
  string snode_size; 
  string sparticle_latspacing;


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

 // Reading the velocity field and longitude and latitude in the velocity file correponds to start_date
  if(LoadVelocityGrid(reference_date)!=0)
    cout << "Error in reading velocities"<< endl;
  

  string name_VTKfile_bathymetry="bathymetry.vtk";

  MakeVTKStructuredGrid2D(vgrid_lon, vgrid_lat, start_depth, bathymetry, nlon, nlat, name_VTKfile_bathymetry );

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////  CONSTRUCTION of NETWORK GRID ///////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(verbose == 1) // Print screen for debugging
    {
      cout <<"GRID CONSTRUCTION PARAMETERS:"<<endl; 
      cout << " Grid step is:"<< node_size<<endl;
      cout << " Low left grid corner (lat,long)=(" << network_ll_lat << "," << network_ll_lon << ")"<<endl;
      cout << " Top right grid corner (lat,long)=(" << network_tr_lat << "," << network_tr_lon << ")"<<endl;
    }

  // opening output files
  ofstream file_output_network;
  string name_output_network;

  // Naming convention (according to parameters) for output files
  char buffer[50];

  sprintf(buffer, "networkdomain%d", domain_network);
  sdomain_network = buffer;

  sprintf(buffer, "_vflow%d", velocity_field);
  svelocity_field=buffer;

  sdepth="_depth";

  snode_size="_nodesize";

  name_output_network=
    output_dir + 
    sdomain_network + 
    svelocity_field + 
    sdepth + DoubletoInt(4, 0, start_depth) +
    snode_size + DoubletoInt(4, 3,node_size) + 
    postfixgrid;
  
  file_output_network.open(name_output_network.c_str());

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
  
  capacity_node = ((network_tr_lat-network_ll_lat)*(network_tr_lon-network_ll_lon)) /(node_size*node_size);

  upnode.reserve(capacity_node);//reserve space for allocation, it reduces the execution time

  // Top right corner of network grid = (network_tr_lat,network_tr_lon)
  // Low left corner of network grid = (network_ll_lat,network_ll_lon)
  // Top right center of network grid = (latmax,lonmax)
  // Low left center of network grid = (latmin,lonmin) 
  delta_lat = node_size;
  latmin = network_ll_lat + (delta_lat/2.0);
  latmax = network_tr_lat - (delta_lat/2.0);
  for(lat=latmin; lat<latmax; lat+=delta_lat)
    {
      delta_lon=node_size/(cos(lat*pig180));
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
	  upnode[num_nodes].center.z = start_depth;

	  // finding two corners of the node from the center
	  upnode[num_nodes].ll.x = lon - delta_lon/2.0;
	  upnode[num_nodes].ll.y = lat - delta_lat/2.0;
	  upnode[num_nodes].ll.z = start_depth;

	  upnode[num_nodes].tr.x = lon + delta_lon/2.0;
	  upnode[num_nodes].tr.y = lat + delta_lat/2.0;
	  upnode[num_nodes].tr.z = start_depth;
	  
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
		  if(start_depth < bathymetry[i][j])
		    { 
		      land++;
		    }
		  point++;
		}
	    }
	  
	  // compute proportion of land per node (land_ratio=1 means land only)
	  upnode[num_nodes].land_ratio = ((double) land)/((double) point);
			    
	  // For each node (1 node = 1 line in output file): write lon, lat of node center + longitude extension + land ratio
	  file_output_network << upnode[num_nodes].center.x << " " << upnode[num_nodes].center.y << " " << delta_lon << " " << upnode[num_nodes].land_ratio << endl;
	  num_nodes++;
	}
    }

  file_output_network.close();

  // Print out total number of nodes if verbose = 1 
  if(verbose == 1)
    {
      cout <<" Total number of nodes: "<< num_nodes << endl;   
      cout <<" Total number of nodes <vector>: "<< upnode.size() << endl;   
      cout <<" Capacity nodes <vector>: "<< upnode.capacity() << endl;   
    }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////  CONSTRUCTION of FINAL NETWORK //////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // opening output files
  ofstream file_output_final_network;
  string name_output_final_network;

  name_output_final_network=
    output_dir + 
    sdomain_network + 
    svelocity_field + 
    sdepth + DoubletoInt(4, 0, final_depth) +
    snode_size + DoubletoInt(4, 3,node_size) + 
    postfixgrid;
  
  file_output_final_network.open(name_output_final_network.c_str());

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
      
      upnode[num_nodes].center.z = final_depth;
      upnode[num_nodes].ll.z = final_depth;
      upnode[num_nodes].tr.z = final_depth;
      
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
		  if(final_depth < bathymetry[i][j])
		    { 
		      land++;
		    }
		  point++;
		}
	    }
	  
	  // compute proportion of land per node (land_ratio=1 means land only)
	  downnode[q].land_ratio = ((double) land)/((double) point);
			    
	  // For each node (1 node = 1 line in output file): write lon, lat of node center + longitude extension + land ratio
	  file_output_final_network << downnode[q].center.x << " " << downnode[q].center.y << " " <<  downnode[q].tr.x-downnode[q].ll.x << " " << downnode[q].land_ratio << endl;
    }

  file_output_final_network.close();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////  CONSTRUCTION of INITIAL POSITIONS OF TRACERS //////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // opening output files
  ofstream file_output_tracer;
  string name_output_tracer;

  sparticle_latspacing = "_particlelatspacing";
 
  //Name output file for initial positions of tracers
  name_output_tracer = 
    output_dir + 
    sdomain_network +
    svelocity_field + 
    sdepth + DoubletoInt(4, 0, start_depth) +
    sparticle_latspacing + DoubletoInt(5,4,particle_latspacing) + 
    postfixitracer;
  
  file_output_tracer.open(name_output_tracer.c_str());

  ////////////////////////////////////////////////////////////////////////
  ////// Computing tracer initial position with 2 loops //////////////////
  ////////////////////////////////////////////////////////////////////////

  // Start/end seeding tracers
  latmin = network_ll_lat;
  latmax = network_tr_lat;
  delta_lat = particle_latspacing;

  lonmin = network_ll_lon;
  lonmax = network_tr_lon;

  if(verbose == 1) // Print screen for debugging
    {
      cout <<"TRACERS CONSTRUCTION PARAMETERS:"<<endl; 
      cout << " Particle latspacing:"<< particle_latspacing <<endl;
      cout << " Low left grid corner (lat,long)=(" << latmin << "," << lonmin << ")"<<endl;
      cout << " Top right grid corner (lat,long)=(" << latmax << "," << lonmax << ")"<<endl;
    }

  int num_tracers = 0;
  for(lat=latmin; lat<latmax; lat+=delta_lat)
    {
      
      delta_lon=particle_latspacing/(cos(lat*pig180));         

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
	  if((start_depth >= bathymetry[i_ll][j_ll]) || 
	     (start_depth >= bathymetry[i_ll][j_tr]) || 
	     (start_depth >= bathymetry[i_tr][j_ll]) || 
	     (start_depth >= bathymetry[i_tr][j_tr]))
	    {
	      // Print out initial positions (lat, lon) of tracers (1 line per particle in output files) 
	      file_output_tracer << lon << " " << lat << endl;
	      num_tracers++;
	    }
	}
    }
  file_output_tracer.close();
  // Print out total number of particles if verbose = 1 
  if(verbose==1)
    {
      cout <<" Total number of tracers: "<< num_tracers << endl;
    } 

  FreeMemoryVelocityGrid();

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

