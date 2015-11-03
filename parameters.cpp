#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>  
#include <cmath>  

using namespace std;
#include "parameters.h"

enum { /* 0 */ DOMAIN_NETWORK,
       /* 1 */ VELOCITY_FIELD,
       /* 2 */ VERBOSE,
       /* 3 */ NODE_SIZE,
       /* 4 */ PARTICLE_LATSPACING,
       /* 5 */ START_DEPTH,
       /* 6 */ FINAL_DEPTH,       
       /* 7 */ START_DATE,
       /* 8 */ INT_STEP,
       /* 9  */ V_SINKING,
       NPARAMETERS
};

int domain_network;
int land_sea_norm;
int velocity_field;
int verbose;
double node_size;
double particle_latspacing;
double start_depth;
double final_depth;
date start_date;
int tau;// it not needed to give a value
double int_step;
double v_sinking;

/* Network  domain:*/

double network_ll_lat;
double network_ll_lon;
double network_tr_lat;
double network_tr_lon;

/* Constants*/
const double pig = 3.14159265358;
const double pig180 = (pig/180.);

/* Configuration parameters*/

const char output_dir[] = "./";
const char input_dir[] = "./";

// File names Outputs
const string postfixgrid = ".grid";
const string postfixitracer = ".itrac";
const string postfixftracer = ".ftrac";
const string postfixmatrix = ".matrix";

char const ncdir[] = "/data/geo/escola/roms_benguela/";
date reference_date = {8,  //year
		       1,  //month
		       1,  //day
};
void listofparameters(string pname[], string ptype[])
{
  int parameter;

  cout << "In the file parameters there must be these parameters ";
  cout << "(it is not necessary the same order, but it is recommended):"<< endl;
  cout << "Number Name Type"<< endl;
  for(parameter=0; parameter<NPARAMETERS; parameter++)
    cout << parameter <<"->"<<pname[parameter]<<" "<<ptype[parameter]<<endl;
}

int readparams(string nfileparameters)
{
  string me="readparams()";
  ifstream fparameters(nfileparameters.c_str());

  int delimiter,end;
  string name, value;
  string line;


  int *pflag;

  int parameter;
  string pname[NPARAMETERS];
  string ptype[NPARAMETERS];
  
  if (!fparameters.is_open())
    {
      cout << me <<": Skipping unreadable file \"" << nfileparameters.c_str() << "\" "<<endl; 
      return 1;
    }

  pflag = (int*) calloc (NPARAMETERS,sizeof(int));

  pname[DOMAIN_NETWORK]="domain_network";           ptype[DOMAIN_NETWORK]="int";
  pname[VELOCITY_FIELD]="velocity_field";           ptype[VELOCITY_FIELD]="int"; 
  pname[VERBOSE]="verbose";                         ptype[VERBOSE]="int"; 
  pname[NODE_SIZE]="node_size";                     ptype[NODE_SIZE]="double"; 
  pname[PARTICLE_LATSPACING]="particle_latspacing"; ptype[PARTICLE_LATSPACING]="double"; 
  pname[START_DEPTH]="start_depth";                 ptype[START_DEPTH]="double"; 
  pname[FINAL_DEPTH]="final_depth";                 ptype[FINAL_DEPTH]="double"; 
  pname[START_DATE]="start_date";                   ptype[START_DATE]="int-int-int";
  pname[INT_STEP]="int_step";                       ptype[INT_STEP]="double"; 
  pname[V_SINKING]="vsinking";                     ptype[V_SINKING]="double"; 

  while(!fparameters.eof())
    {
      getline(fparameters,line);
      if (line[0] == '#') 
	continue;  /* ignore comment line which starts with #*/

      delimiter = line.find("=");
      end = line.length();
      if (end == 0) 
	continue; /* ignore blank line */
      
      value = line.substr(delimiter+1, end);
      name = line.substr(0, delimiter);

      for(parameter=0; parameter<NPARAMETERS; parameter++)
	{
	  if (name.compare(pname[parameter]) == 0)
	    {
	      pflag[parameter]++;
	      if(pflag[parameter]>1)
		{
		  cout  << me << ": Parameter "<< name << " repeated"<<endl;
		  listofparameters(pname, ptype);
		  return 1;
		}
	       if ((delimiter+1) == end) 
		 {
		   cout << me << ": Parameter "<< name << " has no value"<<endl;
		   listofparameters(pname, ptype);
		   return 1;
		 }
	       switch (parameter) 
		{
		case DOMAIN_NETWORK:{
		  domain_network = atoi(value.c_str());}
		  break;
		case VELOCITY_FIELD:
                  velocity_field = atoi(value.c_str());
		  break;
		case VERBOSE:
		  verbose = atoi(value.c_str());
		  break;
		case NODE_SIZE:
		  node_size = atof(value.c_str());
		  break;
		case PARTICLE_LATSPACING:
		  particle_latspacing = atof(value.c_str());
		  break;
		case START_DEPTH:
		  start_depth = atof(value.c_str());
		  start_depth = (-1.0)*start_depth;
		  break;
		case FINAL_DEPTH:
		  final_depth = atof(value.c_str());
		  final_depth = (-1.0)*final_depth;
		  break;
		case START_DATE:
		  if(sscanf(value.c_str(),"%u-%u-%u",&start_date.day,&start_date.month,&start_date.year)!=3)
		    {
		      cout  << me << ": Date format in start_date is incorrect" << endl;
		      listofparameters(pname, ptype);
		      return 1;
		    }
		  break;
		case INT_STEP:
		  int_step = atof(value.c_str());
		  break;
		case V_SINKING:
		  v_sinking = atof(value.c_str());
		  break;  
		default:
		  cout  << me << ": Unknown parameter "<< name <<endl ;
		  listofparameters(pname, ptype);
		  return 1;
		}
	       break;
	    }
	}
      if(parameter==NPARAMETERS)
	{
	  cout  << me << ": Unknown parameter "<< name <<endl ;
	  listofparameters(pname, ptype);
	  return 1;
	}
    }

  for(parameter=0; parameter<NPARAMETERS; parameter++)
	{
	  if(pflag[parameter]==0)
	    {
	      cout  << me << ": parameter "<< pname[parameter] << " is not defined"<<endl;
	      listofparameters(pname, ptype);
	      return 1;
	    }
	}

  tau = 1 + (int) ((2*fabs(final_depth-start_depth))/(v_sinking));
  /* VERBOSE */
  if(verbose == 1)
    {
      cout << "PARAMETERS FROM FILE: "<< nfileparameters <<endl; 
      cout << " domain_network = "<<domain_network << endl; 
      cout << " velocity_field = "<<velocity_field<< endl ;
      cout << " verbose = "<< verbose<< endl;
      cout << " node_size = "<< node_size<< endl;
      cout << " particle_latspacing = "<<particle_latspacing<<endl ;
      cout << " start_depth = "<<start_depth<<endl;
      cout << " final_depth = "<<final_depth<<endl;
      cout << " start_date = "<< start_date.day<<"-"<<start_date.month<<"-"<<start_date.year<<endl ;
      cout << " tau = "<<tau<< endl ;
      cout << " int_step = "<<int_step<< endl ;
      cout << " v_sinking = "<<v_sinking<< endl ;
    }
  
  fparameters.close();

  switch(velocity_field) {
  case 1:
    {
      //degree_resolution = 0.0625; // 1/16 degrees
    }
    break;
    default: 
       {
	 cout << "Incorrect value of velocity_field variable"<<endl ;
	 listofparameters(pname, ptype);
	 return 1;
       }
  }


  /* DOMAIN */

  enum { /* 1 */ WHOLE_DOMAIN=1,
	 /* 2 */ INNER_DOMAIN,
	 /* 3 */ NORTH_SUBDOMAIN,
	 /* 4 */ SOUTH_SUBDOMAIN,
  };

    switch(domain_network) {
    case WHOLE_DOMAIN:
      {
	network_ll_lat= -35.6495;
	network_ll_lon= 3.83333;
	network_tr_lat= -12.0508;
	network_tr_lon= 19.9167;
      }
      break;
    case INNER_DOMAIN:
      {
	network_ll_lat=-35.00;
	network_ll_lon=4.0;
	network_tr_lat=-13.0;
	network_tr_lon=19.00;
      }
      break;
    case NORTH_SUBDOMAIN:
      {
	network_ll_lat=-27.0;
	network_ll_lon=4.0;
	network_tr_lat=-13.0;
	network_tr_lon=19.0;
      }
      break;
    case SOUTH_SUBDOMAIN:
      {
	network_ll_lat=-35.0;
	network_ll_lon=4.0;
	network_tr_lat=-27.0;
	network_tr_lon=19.0;
      }
      break;
     default: 
       {
	 cout << "Incorrect value of domain_network variable"<<endl ;
	 listofparameters(pname, ptype);
	 return 1;
       }
    }

    //deallocate memory:

    free (pflag);
    return 0;
}
