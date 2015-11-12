#include <iostream>
#include <iomanip>
#include <ctime>  
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;
#include "parameters.h" // Function to read parameters.dat
#include "ioutil.h" 
#include "vectorXYZ.h" 

/////////////////////////////////////////////////////////////////////////////////////////////////
// Declaring constants and extern variables
/////////////////////////////////////////////////////////////////////////////////////////////////

extern int verbose;
extern double nodesize;
extern double startdepth;
extern double finaldepth;
extern int tau;
extern double network_ll_lat;
extern double network_ll_lon;
extern double network_tr_lat;
extern double network_tr_lon;

extern string filename_upnetwork;
extern string filename_downnetwork;
extern string filename_itracer;
extern string filename_ftracer;
extern string filename_matrix;

double pig=acos(-1.0);
double pig180=(pig/180.0);

/////////////////////////////////////

struct rectangle {
  struct vectorXYZ center;
  struct vectorXYZ tr;
  struct vectorXYZ ll;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////// Prototype-functions used in the main code (written below, at the end of the main code)
/////////////////////////////////////////////////////////////////////////////////////////////////

void print_usage(const string me);

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
      cout << "Error in reading file parameters" << endl;// Show a error message 
      return 1;
    }

  //////////////////// READ TRACERS INITIAL POSITIONS //////////////////

  // Count number of lines = total number of lagrangian particles
  int num_itracers=CountLines(filename_itracer.c_str());
  if(num_itracers < 0)
    {
      cout << "ERROR in counting lines" << endl;
      return 1;
    }

  // Print initial number of tracers to screen for debugging
  if(verbose == 1) 
    {
      cout <<" Initial number of initial tracers = "<< num_itracers << endl;
    }

  double *itracer_lon;
  double *itracer_lat;

  itracer_lon = new double [num_itracers];
  itracer_lat = new double [num_itracers];

  // open tracer file 
  ifstream ifile_itracer(filename_itracer.c_str());// open tracer file 
  
  int i;

  // Read initial positions of particles 
  for(i=0; i<num_itracers; i++)
    {
      ifile_itracer >> itracer_lon[i] >> itracer_lat[i];
    }
  ifile_itracer.close();

 //////////////////// READ TRACERS FINAL POSITIONS //////////////////
  
  // Count number of lines = total number of lagrangian particles
  int num_ftracers=CountLines(filename_ftracer.c_str());
  if(num_ftracers < 0)
    {
      cout << "ERROR in counting lines" << endl;
      return 1;
    }

  // Test to know if initial and final number of particle is conserved
  if(num_itracers != num_ftracers)
    cout <<" initial number of tracers is not equal to final number "<< endl;

  // If initial and final number agree, keep only 1 variable "num_tracers"
  cout <<" Final number of tracers = "<< num_ftracers << endl;
  
  int num_tracers;

  num_tracers = num_ftracers;

  double *ftracer_lon;
  double *ftracer_lat;
  double *ftracer_depth;
  double *timespent;

  ftracer_lon = new double [num_ftracers];
  ftracer_lat = new double [num_ftracers];
  ftracer_depth = new double [num_ftracers];
  timespent = new double [num_ftracers];

  // open tracer file and read final positions 
  ifstream ifile_ftracer(filename_ftracer.c_str());
  
  for(i=0; i<num_ftracers; i++)
    {
      ifile_ftracer >> ftracer_lon[i] >> ftracer_lat[i] >> ftracer_depth[i]  >> timespent[i] ;
    }
  ifile_ftracer.close();

  //////////////////// READ NETWORK GRID //////////////////

  ifstream ifile_upnetwork(filename_upnetwork.c_str());
  ifstream ifile_downnetwork(filename_downnetwork.c_str());

  // Count number of lines = total number of nodes
  int num_nodes=CountLines(filename_upnetwork.c_str());
  if(num_nodes < 0)
    {
      cout << "ERROR in counting lines" << endl;
      return 1;
    }

  // Print total number of nodes to screen for debugging
  if(verbose == 1)
    {
      cout <<" Total number of nodes = "<< num_nodes << endl;
    }

  rectangle *node;
  double *upland_ratio, *downland_ratio;
  node = new rectangle [num_nodes];
  upland_ratio = new double [num_nodes];
  downland_ratio = new double [num_nodes];


  ///////////////////////////////////// NODE INDICES CALCULATION ////////////////////
  double delta_lon;
  double *land_ratio;
  double delta_node_lat;
  double *delta_node_lon;
  int **node_index;
  int max_iindex_node;
  double max_node_lat;
  int *max_jindex_node;
  double *max_node_lon;
  int count,j;
  
  // Constant latitudinal widht of nodes
  delta_node_lat = nodesize;
  
  // Compute node index (integer) maximum (upper limits) for latitude
  max_iindex_node = (int) ((network_tr_lat-network_ll_lat)/delta_node_lat);
  if(0.0==fmod(network_tr_lat-network_ll_lat-delta_node_lat,delta_node_lat))
     max_iindex_node--;// manual fixing a little error!! 

  // Compute actual latitude (double) corresponding to the upper limit of the last node (can be slightly different than network_tr_lat)
  max_node_lat = delta_node_lat * ((double)max_iindex_node ) + network_ll_lat;
  
  // Print latitudinal properties (max index & max lat) of nodes to screen for debugging
  if(verbose == 1) 
    {
      cout << " Maximum latitudinal index (node) = " << max_iindex_node << endl;
      cout << " Maximum latitudinal limit (top of the node) = " << max_node_lat << endl;
    }
  
  // Declaration of variables (!!Some are irregular matrices!!)
  delta_node_lon =  new double [max_iindex_node];
  max_node_lon =  new double [max_iindex_node];
  max_jindex_node =  new int [max_iindex_node];
  node_index =  new int *[max_iindex_node];
  land_ratio = new double [num_nodes];
  int maxmax_jindex_node = 0;
  double maxmax_node_lon = 0;
  
  // Initialize counter (equivalent to the absolute number of each node)
  count = 0;

  for (i = 0; i < max_iindex_node; i++)
    {
      // Read (line by line) of network grid file
      ifile_upnetwork >> node[count].center.x >> node[count].center.y >> delta_lon >> upland_ratio[count];
      ifile_downnetwork >> node[count].center.x >> node[count].center.y >> delta_lon >> downland_ratio[count];

      node[count].center.z = startdepth;
      
      // finding two corners of the node from the center
      node[count].ll.x = node[count].center.x - delta_lon/2.0;
      node[count].ll.y = node[count].center.y - nodesize/2.0;
      node[count].ll.z = startdepth;
      
      node[count].tr.x = node[count].center.x + delta_lon/2.0;
      node[count].tr.y = node[count].center.y + nodesize/2.0;
      node[count].tr.z = startdepth;
      
      // Latitudinal line per latitudinal line (because irregular grid in longitude)
      delta_node_lon[i] = delta_lon;
      max_jindex_node[i] =  (int) ((network_tr_lon-network_ll_lon)/delta_node_lon[i]);
      if(0.0==fmod(network_tr_lon-network_ll_lon-delta_node_lon[i],delta_node_lon[i]))
	max_jindex_node[i]--;// manual fixing a little error!! 
      
      max_node_lon[i] = delta_node_lon[i]*((double) max_jindex_node[i]) + network_ll_lon;
      
      // Check values and save the maximum of both (for debugging)
      if(max_jindex_node[i] > maxmax_jindex_node)
	maxmax_jindex_node = max_jindex_node[i];

      if(max_node_lon[i] > maxmax_node_lon)
	maxmax_node_lon = max_node_lon[i];
      
      // Allocate memory
      node_index[i] = new int [max_jindex_node[i]];
      // First node of the line at latitude i
      node_index[i][0] = count;
      count++;

      // Loop from 1 to j < max_jindex_node[i] to read all lines that have same delta_lon = same latitude = same i index (and only them)
      for (j = 1; j < max_jindex_node[i]; j++)
	{
	  // Read (another line) of network grid file
	  ifile_upnetwork >> node[count].center.x >> node[count].center.y >> delta_lon >> upland_ratio[count];
	  ifile_downnetwork >> node[count].center.x >> node[count].center.y >> delta_lon >> downland_ratio[count];

	  node[count].center.z = startdepth;
	  
	  // finding two corners of the node from the center
	  node[count].ll.x = node[count].center.x - delta_lon/2.0;
	  node[count].ll.y = node[count].center.y - nodesize/2.0;
	  node[count].ll.z = startdepth;
	  
	  node[count].tr.x = node[count].center.x + delta_lon/2.0;
	  node[count].tr.y = node[count].center.y + nodesize/2.0;
	  node[count].tr.z = startdepth;
	  
	  node_index[i][j] = count;
	  // Add 1 at the end because node index starts at 0
	  count++;
	}
    }

  // Print longitudinal properties (max index & max lon) of nodes to screen for debugging
  if(verbose == 1) 
    {
      cout << " Maximum longitudinal index (node) = " << maxmax_jindex_node << endl;
      cout << " Maximum longitudinal limit (top of the node) = " << maxmax_node_lon << endl;	  
    }

  // Print count (should equal total number of nodes) to screen for debugging
  if(verbose == 1) 
    {
      cout << " Maximum values of node index (should equal total number of nodes) = " << count << endl;
    }

  // Close file
  ifile_upnetwork.close();
  ifile_downnetwork.close();

  ///////////////////////////////////// TRANSPORT MATRIX CALCULATION ////////////////////////////
  cout << "TRANSPORT MATRIX CALCULATION" << endl;
  
  int initial_node, final_node;
  int **trans_matrix;
  double **sum_timespent,**sum2_timespent; 
  double  min_node_lat, min_node_lon;

  // Declaring size of matrix (allocating memory)
  trans_matrix = new int *[num_nodes];
  sum_timespent = new double *[num_nodes];
  sum2_timespent = new double *[num_nodes];
  for ( i=0;i<num_nodes;i++)
    {
      trans_matrix[i] = new int [num_nodes];
      sum_timespent[i] = new double [num_nodes];
      sum2_timespent[i] = new double [num_nodes];
    }

  // Initialize transport matrix to 0  
  for ( i = 0; i < num_nodes; i++)
    for ( j = 0; j < num_nodes; j++)
      {
	trans_matrix[i][j] = 0;
	sum_timespent[i][j] = 0.0;
	sum2_timespent[i][j] = 0.0;
      }
  // Lower left limits of the domain (to check if particle is inside)
  min_node_lat = network_ll_lat;
  min_node_lon = network_ll_lon;

  int num_blind_particles=0;
  int num_outer_lat_particles=0;
  int num_outer_lon_particles=0;
  int num_sticky_particles=0;
  int num_slow_particles=0;

  // Loop over all particles
  for(int trac_index = 0; trac_index < num_tracers; trac_index++)
    { 
      ///////// INITIAL POSITION //////////////////////////////////////
      // Check if the particle is within latitudinal limits    
      if((itracer_lat[trac_index]>= min_node_lat)&&(itracer_lat[trac_index]< max_node_lat))
	// Compute the corresponding latitudinal index of the node where the particle is
	i = (int) ((itracer_lat[trac_index] - min_node_lat)/delta_node_lat);
      else
	{
	  num_blind_particles++;
	  continue;//these are the blind particles
	}
      // Check if the particle is within longitudinal limits    
      if((itracer_lon[trac_index]>= min_node_lon)&&(itracer_lon[trac_index] < max_node_lon[i]))
	// Compute the corresponding longitudinal index of the node where the particle is
	j = (int) ((itracer_lon[trac_index] - min_node_lon)/delta_node_lon[i]);
      else
	{
	  num_blind_particles++;
	  continue;//these are the blind particles
	}
      // Go from 2 indices [i],[j] to absolute number of node (initial position)
      initial_node = node_index[i][j];

      ///////// FINAL POSITION //////////////////////////////////////

      if((timespent[trac_index]>tau))
	{
	  num_slow_particles++;
	  continue;// these are the slow particles
	}


      // Check if the particle reached the final depth 
      if((ftracer_depth[trac_index]>finaldepth))
	{
	  if((ftracer_lat[trac_index]<= network_ll_lat)||(ftracer_lat[trac_index]>=network_tr_lat))
	    {
	      num_outer_lat_particles++;
	      continue;
	    }
	  if((ftracer_lon[trac_index]<= network_ll_lon)||(ftracer_lon[trac_index]>=network_tr_lon)) 
	    {
	      num_outer_lon_particles++;
	      continue;
	    }
	  num_sticky_particles++;
	  continue;// these are the outer particles plus sticky
	}

      // Check if the particle is within latitudinal limits    
      if((ftracer_lat[trac_index]>= min_node_lat)&&(ftracer_lat[trac_index]< max_node_lat))
	// Compute the corresponding latitudinal index of the node where the particle is
	i = (int) ((ftracer_lat[trac_index] - min_node_lat)/delta_node_lat);
      else
	{
	  num_blind_particles++;
	  continue;
	}
      // Check if the particle is within longitudinal limits    
      if((ftracer_lon[trac_index]>= min_node_lon)&&(ftracer_lon[trac_index]< max_node_lon[i]))
	// Compute the corresponding longitudinal index of the node where the particle is
	j = (int) ((ftracer_lon[trac_index]-min_node_lon)/delta_node_lon[i]);
      else
	{
	  num_blind_particles++;
	  continue;
	}
      // Go from 2 indices [i],[j] to absolute number of node (final position)
      final_node = node_index[i][j];
      
      // Add 1 to the corresponding elements of the matrix 
      trans_matrix[initial_node][final_node] += 1;
      sum_timespent[initial_node][final_node] += timespent[trac_index];
      sum2_timespent[initial_node][final_node] += (timespent[trac_index]*timespent[trac_index]);
    }
  
  ////////////////////////////////// OUTPUT WRITING ///////////////////////


  int sum_tracers = 0; 
  // Write only non-null elements of the matrix in the output file

  double average_time,variance_time;
  for(int i=0; i<num_nodes; i++)
    {
      for(int j=0; j<num_nodes; j++)
	{
	  if(trans_matrix[i][j]>0)
	    {
	      average_time = sum_timespent[i][j]/((double) trans_matrix[i][j]);
	      variance_time = (sum2_timespent[i][j]/((double) trans_matrix[i][j])) - (average_time*average_time);
	      sum_tracers+=trans_matrix[i][j];
	    }
	}
    }
  
  // Print to screen sum of all elements of the matrix
  if(verbose == 1) 
    {
      cout << " Total sum of matrix elements (should equal total number of tracers) = " << sum_tracers <<endl;
      cout << " Total sum of blind particles = " << num_blind_particles <<endl;
      cout << " Total sum of slow particles = " << num_slow_particles <<endl;
      cout << " Total sum of outer lat particles = " << num_outer_lat_particles <<endl;
      cout << " Total sum of outer lon particles = " << num_outer_lon_particles <<endl;
      cout << " Total sum of sticky particles = " << num_sticky_particles <<endl;
    }

  //the vtk files

  int released_particles;
  int out_degree;

  string filenamegridvtk = filename_matrix +  "Upnode.vtk";

  ofstream fgridvtk(filenamegridvtk.c_str());// output vtk file 
  fgridvtk << "# vtk DataFile Version 3.0" << endl;
  fgridvtk <<  "vtk output" << endl;
  fgridvtk <<  "ASCII " << endl;
  fgridvtk << "DATASET POLYDATA" << endl;
  fgridvtk << "POINTS "<< 4*num_nodes <<" float" << endl;
  for(i = 0; i<num_nodes; i++ )
    {
      fgridvtk << node[i].ll.x <<" "<< node[i].ll.y <<" "<< startdepth*0.001 << endl;
      fgridvtk << node[i].ll.x <<" "<< node[i].tr.y <<" "<< startdepth*0.001 << endl;
      fgridvtk << node[i].tr.x <<" "<< node[i].tr.y <<" "<< startdepth*0.001 << endl;
      fgridvtk << node[i].tr.x <<" "<< node[i].ll.y <<" "<< startdepth*0.001 << endl;
    }
  
  fgridvtk << "POLYGONS  "<< num_nodes <<" "<< 5*num_nodes <<endl;
  for(i = 0; i<num_nodes; i++)
    {
      fgridvtk <<"4 "<< 4*i <<" "<< 4*i+1<<" "<< 4*i+2<<" "<<4*i+3<< endl;
    }
  fgridvtk << "CELL_DATA "<< num_nodes <<endl;
  fgridvtk << "SCALARS upland_ratio float"<<endl;
  fgridvtk << "LOOKUP_TABLE default"<< endl;
  for(i = 0; i<num_nodes; i++)
    {
     fgridvtk << upland_ratio[i] << endl;
    }
  fgridvtk << "SCALARS released_particles int"<<endl;
  fgridvtk << "LOOKUP_TABLE default"<< endl;
  for(i = 0; i<num_nodes; i++)
    {
      released_particles=0;
      for(j=0; j<num_nodes; j++)
	{
	  if(trans_matrix[i][j]>0)
	    {
	      released_particles+=trans_matrix[i][j];	
	    }	
	}
      fgridvtk << released_particles << endl;
    }
  fgridvtk << "SCALARS out_degree int"<<endl;
  fgridvtk << "LOOKUP_TABLE default"<< endl;
  for(i = 0; i<num_nodes; i++)
    {
      out_degree=0;
      for(j=0; j<num_nodes; j++)
	{
	  if(i!=j && trans_matrix[i][j]>0)
	    {
	      out_degree+=1;	
	    }	
	}
      fgridvtk << out_degree << endl;
    }	
  fgridvtk.close();

  int settled_particles;
  int in_degree;

  string filenamedowngridvtk = filename_matrix +  "Downnode.vtk";

  ofstream fdowngridvtk(filenamedowngridvtk.c_str());// output vtk file 

  fdowngridvtk << "# vtk DataFile Version 3.0" << endl;
  fdowngridvtk <<  "vtk output" << endl;
  fdowngridvtk <<  "ASCII " << endl;
  fdowngridvtk << "DATASET POLYDATA" << endl;
  fdowngridvtk << "POINTS "<< 4*num_nodes <<" float" << endl;
  for(i = 0; i<num_nodes; i++ )
    {
      fdowngridvtk << node[i].ll.x <<" "<< node[i].ll.y <<" "<< finaldepth*0.001 << endl;
      fdowngridvtk << node[i].ll.x <<" "<< node[i].tr.y <<" "<< finaldepth*0.001 << endl;
      fdowngridvtk << node[i].tr.x <<" "<< node[i].tr.y <<" "<< finaldepth*0.001 << endl;
      fdowngridvtk << node[i].tr.x <<" "<< node[i].ll.y <<" "<< finaldepth*0.001 << endl;
    }
  
  fdowngridvtk << "POLYGONS  "<< num_nodes <<" "<< 5*num_nodes <<endl;
  for(i = 0; i<num_nodes; i++)
    {
      fdowngridvtk <<"4 "<< 4*i <<" "<< 4*i+1<<" "<< 4*i+2<<" "<<4*i+3<< endl;
    }
  fdowngridvtk << "CELL_DATA "<< num_nodes <<endl;
  fdowngridvtk << "SCALARS downland_ratio float"<<endl;
  fdowngridvtk << "LOOKUP_TABLE default"<< endl;
  for(i = 0; i<num_nodes; i++)
    {
     fdowngridvtk << downland_ratio[i] << endl;
    }
  fdowngridvtk << "SCALARS settled_particles int"<<endl;
  fdowngridvtk << "LOOKUP_TABLE default"<< endl;
  for(i = 0; i<num_nodes; i++)
    {
      settled_particles=0.0;
      for(j=0; j<num_nodes; j++)
	{
	  if(trans_matrix[j][i]>0)
	    {
	      settled_particles+=trans_matrix[j][i];
	    }	
	}
      fdowngridvtk << settled_particles << endl;
      
    }
  fdowngridvtk << "SCALARS in_degree int"<<endl;
  fdowngridvtk << "LOOKUP_TABLE default"<< endl;
  for(i = 0; i<num_nodes; i++)
    {
      in_degree=0;
      for(j=0; j<num_nodes; j++)
	{
	  if(i!=j && trans_matrix[j][i]>0)
	    {
	      in_degree+=1;	
	    }	
	}
      fdowngridvtk << in_degree << endl;
    }	
  fdowngridvtk.close();

  //DEALLOCATE MEMORY
  delete[] node;
  delete[] upland_ratio;
  delete[] downland_ratio;
  delete[] itracer_lat;
  delete[] itracer_lon;
  delete[] ftracer_lat;
  delete[] ftracer_lon;
  delete[] ftracer_depth;
  delete[] timespent;
  delete[] land_ratio;
  delete[] delta_node_lon;
  delete[] max_node_lon;
  delete[] max_jindex_node;

  for (i=0; i<max_iindex_node; i++)
    {
      delete[] node_index[i];
    }
  delete[] node_index;

  for (int i=0; i<num_nodes; i++)
    {
      delete[] trans_matrix[i];
      delete[] sum_timespent[i];
      delete[] sum2_timespent[i];
    }
  delete[] trans_matrix;
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////// functions used in the main code (declared at the beginning of the main code)
/////////////////////////////////////////////////////////////////////////////////////////////////

// This function prints a little command line manual if there are invalid parameters 
void print_usage(const string me) 
{
  cout << "Usage: " << me << " <file_parameters>" << endl << endl;
}
