#include <iostream>
#include <iomanip>
#include <ctime>  
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
#include "parameters.h"
#include "ioutil.h"

string pname[NPARAMETERS];
string ptype[NPARAMETERS];

int networkdomain;
int vflow;
int verbose;
double nodesize;
double particlespacing;
double startdepth;
double finaldepth;
date startdate;
int tau;// it not needed to give a value
double intstep;
double vsink;

/* Network  domain:*/
double network_ll_lat;
double network_ll_lon;
double network_tr_lat;
double network_tr_lon;

/* Configuration parameters*/
char workingdir[] = "/home/pmonroy/ESCOLA/RUNS/SLFNetwork/";

// File names In/Outputs
string filename_upnetwork;
string filename_downnetwork;
string filename_itracer;
string filename_ftracer;
string filename_matrix;

string infix = "_";
string postfixgrid = ".grid";
string postfixitracer = ".itrac";
string postfixftracer = ".ftrac";
string postfixmatrix = ".matrix";

//char velocitydir[] = "/data/geo/escola/roms_benguela/";
char velocitydir[] = "/scratch/pmonroy/";

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
  
  if (!fparameters.is_open())
    {
      cout << me <<": Skipping unreadable file \"" << nfileparameters.c_str() << "\" "<<endl; 
      return 1;
    }

  pflag = (int*) calloc (NPARAMETERS,sizeof(int));

  pname[NETWORKDOMAIN]="networkdomain";           ptype[NETWORKDOMAIN]="int";
  pname[VFLOW]="vflow";           ptype[VFLOW]="int"; 
  pname[VERBOSE]="verbose";                         ptype[VERBOSE]="int"; 
  pname[NODESIZE]="nodesize";                     ptype[NODESIZE]="double"; 
  pname[PARTICLESPACING]="particlespacing"; ptype[PARTICLESPACING]="double"; 
  pname[STARTDEPTH]="startdepth";                 ptype[STARTDEPTH]="double"; 
  pname[FINALDEPTH]="finaldepth";                 ptype[FINALDEPTH]="double"; 
  pname[STARTDATE]="startdate";                   ptype[STARTDATE]="int-int-int";
  pname[INTSTEP]="intstep";                       ptype[INTSTEP]="double"; 
  pname[VSINK]="vsink";                     ptype[VSINK]="double"; 

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
		case NETWORKDOMAIN:{
		  networkdomain = atoi(value.c_str());}
		  break;
		case VFLOW:
                  vflow = atoi(value.c_str());
		  break;
		case VERBOSE:
		  verbose = atoi(value.c_str());
		  break;
		case NODESIZE:
		  nodesize = atof(value.c_str());
		  break;
		case PARTICLESPACING:
		  particlespacing = atof(value.c_str());
		  break;
		case STARTDEPTH:
		  startdepth = atof(value.c_str());
		  startdepth = (-1.0)*startdepth;
		  break;
		case FINALDEPTH:
		  finaldepth = atof(value.c_str());
		  finaldepth = (-1.0)*finaldepth;
		  break;
		case STARTDATE:
		  if(sscanf(value.c_str(),"%u-%u-%u",&startdate.day,&startdate.month,&startdate.year)!=3)
		    {
		      cout  << me << ": Date format in startdate is incorrect" << endl;
		      listofparameters(pname, ptype);
		      return 1;
		    }
		  break;
		case INTSTEP:
		  intstep = atof(value.c_str());
		  break;
		case VSINK:
		  vsink = atof(value.c_str());
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
  fparameters.close();


  // COMPUTE DERIVED PARAMETERS

  //TAU
  tau = 1 + (int) ((2*fabs(finaldepth-startdepth))/(vsink));

  //NAME FILES

  //Name file for up network 
  string sdepth = pname[STARTDEPTH];
  sdepth.erase(0,5);

  filename_upnetwork=
    workingdir + 
    pname[NETWORKDOMAIN] + DoubletoString(1, 0, (double) networkdomain) + infix +
    pname[VFLOW] + DoubletoString(1,0, (double) vflow) + infix +
    sdepth + DoubletoString(4, 0, startdepth) +  infix +
    pname[NODESIZE] + DoubletoString(4, 3,nodesize) + 
    postfixgrid;

  //Name file for down network 

  filename_downnetwork=
    workingdir + 
    pname[NETWORKDOMAIN] + DoubletoString(1, 0, (double) networkdomain) + infix +
    pname[VFLOW] + DoubletoString(1,0, (double) vflow) + infix +
    sdepth + DoubletoString(4, 0, finaldepth) +  infix +
    pname[NODESIZE] + DoubletoString(4, 3,nodesize) + 
    postfixgrid;

  //Name file for initial positions of tracers

  filename_itracer = 
    workingdir + 
    pname[NETWORKDOMAIN] + DoubletoString(1, 0, (double) networkdomain) + infix +
    pname[VFLOW] +  DoubletoString(1, 0, (double)vflow) + infix +
    sdepth + DoubletoString(4, 0, startdepth) + infix +
    pname[PARTICLESPACING] + DoubletoString(5,4,particlespacing) + 
    postfixitracer;

//Name file for initial positions of tracers

  filename_ftracer = 
    workingdir + 
    pname[NETWORKDOMAIN] + DoubletoString(1,0,(double) networkdomain) + infix +
    pname[VFLOW] +  DoubletoString(1,0,(double) vflow) + infix +
    pname[STARTDEPTH] + DoubletoString(4, 0, startdepth) + infix +
    pname[FINALDEPTH] + DoubletoString(4, 0, finaldepth) + infix +
    pname[PARTICLESPACING] + DoubletoString(5,4,particlespacing) + infix +
    pname[STARTDATE] + DoubletoString(2,0,(double) startdate.day)
                     + DoubletoString(2,0,(double) startdate.month)
                     + DoubletoString(4,0,(double) startdate.year) + infix +
    pname[VSINK]+ DoubletoString(3,0,vsink) + infix +
    pname[INTSTEP] +  DoubletoString(4,2,intstep) +
    postfixftracer;

//Name file for initial positions of tracers

  filename_matrix = 
    workingdir + 
    pname[NETWORKDOMAIN] + DoubletoString(1,0,(double) networkdomain) + infix +
    pname[VFLOW] +  DoubletoString(1,0,(double) vflow) + infix +
    pname[STARTDEPTH] + DoubletoString(4, 0, startdepth) + infix +
    pname[FINALDEPTH] + DoubletoString(4, 0, finaldepth) + infix +
    pname[PARTICLESPACING] + DoubletoString(5,4,particlespacing) + infix +
    pname[NODESIZE] + DoubletoString(4, 3,nodesize) + infix +
    pname[STARTDATE] + DoubletoString(2,0,(double) startdate.day)
                     + DoubletoString(2,0,(double) startdate.month)
                     + DoubletoString(4,0,(double) startdate.year) + infix +
    pname[VSINK]+ DoubletoString(3,0,vsink) + infix +
    pname[INTSTEP] +  DoubletoString(4,2,intstep) +
    postfixmatrix;


  switch(vflow) {
  case 1:
    {
      //degree_resolution = 0.0625; // 1/16 degrees
    }
    break;
    default: 
       {
	 cout << "Incorrect value of vflow variable"<<endl ;
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

    switch(networkdomain) {
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
	 cout << "Incorrect value of networkdomain variable"<<endl ;
	 listofparameters(pname, ptype);
	 return 1;
       }
    }


    /* VERBOSE */
    if(verbose == 1)
      {
	cout << "PARAMETERS FROM FILE: "<< nfileparameters <<endl; 
	cout << " networkdomain = "<<networkdomain << endl; 
	cout << " vflow = "<<vflow<< endl ;
	cout << " verbose = "<< verbose<< endl;
	cout << " nodesize = "<< nodesize<< endl;
	cout << " particlespacing = "<<particlespacing<<endl ;
	cout << " startdepth = "<<startdepth<<endl;
	cout << " finaldepth = "<<finaldepth<<endl;
	cout << " startdate = "<< startdate.day<<"-"<<startdate.month<<"-"<<startdate.year<<endl ;
	cout << " tau = "<<tau<< endl ;
	cout << " intstep = "<<intstep<< endl ;
	cout << " vsink = "<<vsink<< endl ;
      }
    
    //deallocate memory:
    free (pflag);
    return 0;
}
