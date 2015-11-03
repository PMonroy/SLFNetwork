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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////// Prototype-functions used in the main code (written below, at the end of the main code)
/////////////////////////////////////////////////////////////////////////////////////////////////

void print_usage(const string me);
int countlines(string namefile);

string DoubletoInt(int ndigits, int ndecimals, double number);

/////////////////////////////////////////////////////////////////////////////////////////////////
// Declaring constant and global variables
/////////////////////////////////////////////////////////////////////////////////////////////////

//const double pig = 3.14159265358; 
//const double deg2rad = pig / 180.;

//const string input_dir = "/home/merluza/shared/";
//const string output_dir = "/home/merluza/shared/";

// File names Outputs
//const string postfixgrid = ".grid";
//const string postfixtracer = ".trac";
//const string postfixmatrix = ".matrix";

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

  string sdomain_network;
  string svelocity_field;
  string sdepth;
  string sparticle_latspacing;
  string sfile_initial_tracers;
  char buffer[50];

  sprintf(buffer, "networkdomain%d", domain_network);
  sdomain_network = buffer;

  sprintf(buffer, "_vflow%d", velocity_field);
  svelocity_field=buffer;

  sdepth="_depth";

  sparticle_latspacing = "_particlelatspacing";

  // Reconstruct name file
  sfile_initial_tracers = 
    input_dir + 
    sdomain_network +
    svelocity_field + 
    sdepth + DoubletoInt(4,0,start_depth) + 
    sparticle_latspacing + DoubletoInt(5,4,particle_latspacing) +
    postfixitracer; 

  // Count number of lines = total number of lagrangian particles
  int num_itracers=countlines(sfile_initial_tracers);
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
  ifstream file_initial_tracers(sfile_initial_tracers.c_str());// open tracer file 

  int i;

  // Read initial positions of particles 
  for(i=0; i<num_itracers; i++)
    {
      file_initial_tracers >> itracer_lon[i] >> itracer_lat[i];
    }
  file_initial_tracers.close();

 //////////////////// READ TRACERS FINAL POSITIONS //////////////////
  string sstart_date;
  string stau;
  string sint_step;
  string sfile_final_tracers; 

  string sidepth = "_idepth";
  string sfdepth = "_fdepth";
  string svsinking = "_vsinking";

  sstart_date="_startdate";
  sprintf(buffer, "_tau%d", (int) tau);
  stau=buffer;
  sint_step = "_intstep";
  
  // Reconstruct name file
  sfile_final_tracers = 
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
  
  // Count number of lines = total number of lagrangian particles
  int num_ftracers=countlines(sfile_final_tracers);
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
  ifstream file_final_tracers(sfile_final_tracers.c_str());
  
  for(i=0; i<num_ftracers; i++)
    {
      file_final_tracers >> ftracer_lon[i] >> ftracer_lat[i] >> ftracer_depth[i]  >> timespent[i] ;
    }
  file_final_tracers.close();

  //////////////////// READ NETWORK GRID //////////////////
  string sfile_network;
  string snode_size;
  snode_size="_nodesize";
  
  // Reconstruct name file
  sfile_network= 
    output_dir + 
    sdomain_network + 
    svelocity_field + 
    sdepth + DoubletoInt(4, 0, start_depth) +
    snode_size + DoubletoInt(4, 3,node_size) + 
    postfixgrid;

  ifstream file_network(sfile_network.c_str());

  // Count number of lines = total number of nodes
  int num_nodes=countlines(sfile_network);
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

  ///////////////////////////////////// NODE INDICES CALCULATION ////////////////////
  double node_lon;
  double node_lat;
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
  delta_node_lat = node_size;
  
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
      file_network >> node_lon >> node_lat >> delta_lon >> land_ratio[count];

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
	  file_network >> node_lat >> node_lon >> delta_lon >> land_ratio[count];
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
  file_network.close();

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
      if((ftracer_depth[trac_index]>final_depth))
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

  string sfile_matrix;
  
  // Construct name file
  sfile_matrix=
    output_dir + 
    sdomain_network + 
    svelocity_field + 
    sidepth + DoubletoInt(4, 0, start_depth) + 
    sfdepth + DoubletoInt(4, 0, final_depth) + 
    sparticle_latspacing + DoubletoInt(5,4,particle_latspacing) +
    snode_size + DoubletoInt(4, 3,node_size) + 
    sparticle_latspacing + DoubletoInt(5,4,particle_latspacing) +
    sstart_date + DoubletoInt(2,0,start_date.day) + DoubletoInt(2,0,start_date.month) + DoubletoInt(4,0,start_date.year) +
    svsinking + DoubletoInt(3,0,v_sinking) +
    sint_step +  DoubletoInt(4,2,int_step) +
    postfixmatrix;

  ofstream file_matrix(sfile_matrix.c_str());
  int sum_tracers = 0;
 
  // Write only non-null elements of the matrix in the output file
  for(int i=0; i<num_nodes; i++)
    {
      for(int j=0; j<num_nodes; j++)
	{
	  if(trans_matrix[i][j]>0)
	    {
	      file_matrix << i <<" "<< j <<" "<< trans_matrix[i][j]<<" "<<sum_timespent[i][j] <<" "<< sum2_timespent[i][j]<<endl;
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

  // Close file
  file_matrix.close();

  //DEALLOCATE MEMORY
  delete[] itracer_lat;
  delete[] itracer_lon;
  delete[] ftracer_lat;
  delete[] ftracer_lon;
  delete[] ftracer_depth;
  delete[] timespent;
  delete[] land_ratio;
  delete[] delta_node_lon;
  delete[]  max_node_lon;
  delete[]  max_jindex_node;

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
      return -1;
    }
  for ( numlines = 0; getline(ifile, line); ++numlines)
    ;

  ifile.close();
  return numlines;
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

