#include <iostream>
#include <iomanip>
#include <netcdfcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <omp.h>

using namespace std;

#include "velocity.h"
#include "constants.h"

// Return this code to the OS in case of failure.
static const int NC_ERR = 2;

//global variables

double rad_resolution, degree_resolution;

int nlon, nlat, ndepth, ntime;

double *vgrid_lon, *vgrid_lat, ****vgrid_depth;
vectorXYZ ****vfield;
int **land_mask;
double **bathymetry;

int LoadVelocityGrid(date rdate, char velocitydir[])
{
  int i,j,q; // loop indices

  char ncfile[256];
  NcError err(NcError::verbose_nonfatal);

  // Open the first Netcdf file
  sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",velocitydir, rdate.year,rdate.month);
  cout << "Reading the netcdf file: " << ncfile;
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);  
  // Check to see if the file was opened.
  if(!dataFile.is_valid())
    return NC_ERR;
  cout << "->OK" << endl;

  NcDim *ncdim_xi_rho;
  if (!(ncdim_xi_rho = dataFile.get_dim("xi_rho")))
    return NC_ERR;
  nlon = ncdim_xi_rho->size();

  NcDim *ncdim_eta_rho;
  if (!(ncdim_eta_rho = dataFile.get_dim("eta_rho")))
    return NC_ERR;
  nlat = ncdim_eta_rho->size();

  NcDim *ncdim_s_rho;
  if (!(ncdim_s_rho = dataFile.get_dim("s_rho")))
    return NC_ERR;
  ndepth = ncdim_s_rho->size();

  // Degree resolution
  degree_resolution = 0.083333333f;
  rad_resolution = rads*degree_resolution;

  // Get pointers to the latitude and longitude variables. 
  NcVar *lonVar, *latVar;
  if (!(lonVar = dataFile.get_var("lon_rho")))
    return NC_ERR;
  if (!(latVar = dataFile.get_var("lat_rho")))
    return NC_ERR;

  if (!lonVar->set_cur(0,0))
    return NC_ERR;
  if (!latVar->set_cur(0,0))
    return NC_ERR;
    
  vgrid_lon = new double [nlon];
  vgrid_lat = new double [nlat];

  // Get the lat/lon data from the file.
  if (!lonVar->get(vgrid_lon, 1, nlon))
    return NC_ERR;
  if (!latVar->get(vgrid_lat, nlat , 1))
    return NC_ERR;


  // Get pointers to the depth variable.
  NcVar *depthVar;
  if (!(depthVar = dataFile.get_var("depth")))
    return NC_ERR;

  if (!depthVar->set_cur(rdate.day, 0, 0, 0))
    return NC_ERR; 

  double *depth;
  depth = new double [1*1*nlat*nlon];
  if (!depthVar->get(&depth[0], 1, 1, nlat, nlon))
    return NC_ERR;

  // Get pointers to the depth variable.
  NcVar *wVar;
  if (!(wVar = dataFile.get_var("w")))
    return NC_ERR;

  if (!wVar->set_cur(rdate.day, 0, 0, 0))
    return NC_ERR; 

  double *w;
  w = new double [1*1*nlat*nlon];
  if (!wVar->get(&w[0], 1, 1, nlat, nlon))
    return NC_ERR;
  dataFile.close();// close Netcdf file

  // Get the bathymetry and land_mask
  bathymetry = new double *[nlon];
  land_mask = new int *[nlon];
  for (i = 0; i < nlon; i++)
    {
      bathymetry[i] = new double [nlat];
      land_mask[i] = new int [nlat];      
    }

  for(i=0; i<nlon; i++)
    {
      for(j=0; j<nlat; j++)
	{
	  q=j*nlon+ i;
	  if(w[q]==0.0)
	    {
	      land_mask[i][j] = 0;
	      bathymetry[i][j] = 0.0;
	    }	  
	  else
	    {
	      land_mask[i][j] = 1; 
	      bathymetry[i][j] = depth[q];
	    }
	}
    }

  delete[] depth;
  delete[] w;

  return 0;
}

void FreeMemoryVelocityGrid() 
{
  
  delete[] vgrid_lon;
  delete[] vgrid_lat;
  int i;
  for (i = 0; i < nlon; i++)
    {
       delete[] bathymetry[i];
       delete[] land_mask[i];
    }
}

int LoadVelocities(date start_date, int tau, char velocitydir[])
{

  int t,i,j,k;
  
  int nlonu, nlatv;

  char ncfile[256];

  NcError err(NcError::verbose_nonfatal);

  // Open the first Netcdf file
  sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",velocitydir, start_date.year,start_date.month);
  cout << "Reading first netcdf file: " << ncfile;
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);  
  // Check to see if the file was opened.
  if(!dataFile.is_valid())
    return NC_ERR;
  cout << " OK" << endl;

  // Get dimensions of varibles

  NcDim *ncdim_xi_u;
  if (!(ncdim_xi_u = dataFile.get_dim("xi_u")))
    return NC_ERR;
  nlonu = ncdim_xi_u->size();
 
  NcDim *ncdim_eta_v;
  if (!(ncdim_eta_v = dataFile.get_dim("eta_v")))
    return NC_ERR;
  nlatv = ncdim_eta_v->size();

  NcDim *ncdim_time;
  if (!(ncdim_time = dataFile.get_dim("time")))
    return NC_ERR;
  ntime = ncdim_time->size();

  dataFile.close();// close Netcdf file

  // Read the variables that depend on time

  double *h, *u, *v, *w;

  h = new double [(tau+2)*nlon*nlat*ndepth];
  u = new double [(tau+2)*nlonu*nlat*ndepth];
  v = new double [(tau+2)*nlon*nlatv*ndepth];
  w = new double [(tau+2)*nlon*nlat*ndepth];

  int npoints_center,npoints_layer_center;
  int npoints_u,npoints_layer_u;
  int npoints_v,npoints_layer_v;

  npoints_center= nlon*nlat*ndepth;
  npoints_layer_center = nlon*nlat;

  npoints_u= nlonu*nlat*ndepth;
  npoints_layer_u = nlonu*nlat;

  npoints_v= nlon*nlatv*ndepth;
  npoints_layer_v = nlon*nlatv;

  NcVar *depth_var, *u_var, *v_var, *w_var;
  unsigned int start_time, time, final_time;
  date tdate;
  unsigned int count, sum_count;
  
  start_time = DATE_TO_TIME(start_date);
  final_time = start_time + tau  + 2;

  time = start_time;
  sum_count = 0;

  while(time < final_time)
    {
      TIME_TO_DATE(tdate,time);
      
      count = ntime - (tdate.day-1);

      if((time + count) > final_time)
	count = final_time - time;
      
      sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",velocitydir, tdate.year,tdate.month);
      cout << "reading netcdf file: " <<ncfile <<" start day="<< tdate.day <<" number of days=" << count;

      NcFile dataFile(ncfile, NcFile::ReadOnly);

      // Check to see if the file was opened.
      if(!dataFile.is_valid())
	return NC_ERR;
 
      //Read depth
      if (!(depth_var = dataFile.get_var("depth")))
	return NC_ERR;
      if (!depth_var->set_cur(tdate.day-1, 0, 0, 0))
	return NC_ERR; 
      if (!depth_var->get(&h[sum_count*npoints_center], count, ndepth, nlat, nlon))
	return NC_ERR;

      //Read depth u
      if (!(u_var = dataFile.get_var("u")))
	return NC_ERR;
      if (!u_var->set_cur(tdate.day-1, 0, 0, 0))
	return NC_ERR; 
      if (!u_var->get(&u[sum_count*npoints_u], count, ndepth, nlat, nlonu))
	return NC_ERR;

      //Read depth v
      if (!(v_var = dataFile.get_var("v")))
	return NC_ERR;
      if (!v_var->set_cur(tdate.day-1, 0, 0, 0))
	return NC_ERR; 
      if (!v_var->get(&v[sum_count*npoints_v], count, ndepth, nlatv, nlon))
	return NC_ERR;

      //Read w
      if (!(w_var = dataFile.get_var("w")))
	return NC_ERR;
      if (!w_var->set_cur(tdate.day-1, 0, 0, 0))
	return NC_ERR; 
      if (!w_var->get(&w[sum_count*npoints_center], count, ndepth, nlat, nlon))
	return NC_ERR;

      dataFile.close();
      cout << " OK" << endl;
      
      time += count;
      sum_count += count;
    }

  /* Transpose the arrays */

  // Dynamically allocate velocity fields


  vgrid_depth = new double ***[tau+2];
  vfield = new vectorXYZ ***[tau+2];

  for (t = 0; t < (tau+2); t++)
    {
      vgrid_depth[t] = new double **[nlon];
      vfield[t] = new vectorXYZ **[nlon];
      for (i = 0; i < nlon; i++)
	{
	  vgrid_depth[t][i] = new double *[nlat];
	  vfield[t][i] = new vectorXYZ *[nlat];
	  for(j = 0; j < nlat; j++)
	    {
	      vgrid_depth[t][i][j] = new double [ndepth];
	      vfield[t][i][j] = new vectorXYZ [ndepth];
	    }
	}
    }

#pragma omp parallel for default(shared) private(i,j,k)

  for(t = 0; t < tau+2; t++)
    {
      for(i = 0; i < nlon; i++)
	{
	  for(j = 0; j < nlat; j++)
	    {
	      for(k = 0; k < ndepth; k++)
		{
		  vgrid_depth[t][i][j][k] = *(h + t*npoints_center+k*npoints_layer_center+ j*nlon+ i);
		  vfield[t][i][j][k].z= *(w + t*npoints_center+k*npoints_layer_center+ j*nlon+ i);
		  vfield[t][i][j][k].x= (*(u + t*npoints_u+k*npoints_layer_u+ j*nlonu+ i-(i==(nlon-1)))+*(u + t*npoints_u+k*npoints_layer_u+ j*nlonu + i-(i!=0)))/2.0;
		  vfield[t][i][j][k].y= (*(v + t*npoints_v+k*npoints_layer_v+ (j-(j==(nlat-1)))*nlon+i)+*(v + t*npoints_v+k*npoints_layer_v+ (j-(j!=0))*nlon + i))/2.0;

		  vfield[t][i][j][k] = secondsday *  vfield[t][i][j][k];

		}
	    }
	}
    }

  delete[] h;
  delete[] u;
  delete[] v;
  delete[] w;

  return 0;
}

void FreeMemoryVelocities(int tau) 
{
  
  int t,i,j;
   for(t = 0; t < tau+2; t++)
    {
      for(i = 0; i < nlon; i++)
	{
	  for(j = 0; j < nlat; j++)
	    {
	      delete[] vgrid_depth[t][i][j];
	      delete[] vfield[t][i][j];
	    }
	}
    }
}

void locate(double xx[], unsigned long n, double x, unsigned long *j)
{

  /* Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j]
   * and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned
   * to indicate that x is out of range.
   */

  unsigned long ju,jm,jl;
  int ascnd;
  
  jl=0;
  ju=n+1;
  ascnd=(xx[n] > xx[1]);
  while ((ju-jl) > 1) 
    {
      jm=(ju+jl) >> 1;
      if (x > xx[jm] == ascnd)
	jl=jm;
      else
	ju=jm;
    }
  if (x == xx[1]) 
    *j=1;
  else if(x == xx[n]) 
    *j=n-1;
  else 
    *j=jl;
}

int GetLonIndex(double longitude)
{
  int index;

  /* Locate index longitude*/
  index = (int) ((longitude-vgrid_lon[0])/(degree_resolution));

  return index;
}

int GetLatIndex(double latitude)
{
  int index;
  double *llat;
  unsigned long index_latitude;

  /* Locate index latitude*/
  llat = vgrid_lat - 1;
  locate(llat, nlat, latitude, &index_latitude);
  if(index_latitude == 0 || index_latitude == abs(nlat))
    return 1;
  else
    index = index_latitude - 1;

  return index;
}

int GetIndices(unsigned long time, vectorXYZ point, vectorIJK *index)
{

  /* point.x -> longitude
   * point.y -> latitude
   * point.z -> depth
   */

  double *llat,*ddpt ;
  unsigned long index_latitude=0, index_depth=0;
  int i,j;

  /* Locate index longitude*/

  index->i = (int) ((point.x-vgrid_lon[0])/(degree_resolution));

  if(index->i < 0 || index->i >= nlon-1)
      return 1;

  /* Locate index latitude*/
  llat = vgrid_lat - 1;
  locate(llat, nlat, point.y, &index_latitude);
  if(index_latitude == 0 || index_latitude == abs(nlat))
    return 1;
  else
    index->j = index_latitude - 1;
  
  /* Locate index depth */
  for(i = 0; i < 2; i++)
    {
      for(j = 0; j < 2; j++)
	{
	  ddpt = vgrid_depth[time][index->i+i][index->j+j] - 1;
	  locate(ddpt, ndepth, point.z, &index_depth);	     
	  if(index_depth == 0 || index_depth == abs(ndepth))
	    return 1;
	  else
	    index->k[i][j] = index_depth - 1;
	}
    }
  return 0;
}

int GetVelocity(double t,vectorXYZ point, vectorXYZ *vint)
{
  vectorXYZ vgrid[16];
  vectorXYZ vcomp[16];
  vectorIJK index[2];
  unsigned long time;

  double alpha, beta;

  int h,i,j,k;
  int deltatime,deltai,deltaj,deltak;
  int q;

  /* Calculates the vectors vb[15] and points ptmb[15] */
  time = (unsigned long) t;

  if(GetIndices(time, point, &index[0])==1)
    return 1;
  
  if(GetIndices(time+1, point, &index[1])==1)
    return 1;

  /* Vectors and points with only one int index*/
  q=0;  
  for(deltatime=0; deltatime<2; deltatime++)
    {
      for(deltai=0; deltai<2; deltai++)
	{
	  for(deltaj = 0; deltaj<2; deltaj++)
	    {
	      for(deltak=0; deltak<2; deltak++)
		{
		  i = index[deltatime].i + deltai;
		  j = index[deltatime].j + deltaj;
		  k = index[deltatime].k[deltai][deltaj] + deltak;
		  h = time+deltatime;

		  vgrid[q].x = rads*vgrid_lon[i];
		  vgrid[q].y = log(fabs((1.0/cos(rads*vgrid_lat[j]))+tan(rads*vgrid_lat[j])));
		  vgrid[q].z = vgrid_depth[h][i][j][k];
		  
		  vcomp[q] = vfield[h][i][j][k];
		  
		  q++; 
		}
	    }
	}
    }

  point.x = rads*point.x;
  point.y = log(fabs((1.0/cos(rads*point.y))+tan(rads*point.y)));

  /* COAST CHECKING */
  if( (point.z < bathymetry[index[0].i  ][index[0].j  ]) &&
      (point.z < bathymetry[index[0].i+1][index[0].j  ]) &&
      (point.z < bathymetry[index[0].i  ][index[0].j+1]) &&
      (point.z < bathymetry[index[0].i+1][index[0].j+1]))
    return 1;// The particle has reached the coast

  /* Depth Interpolation: */
  for(q=0; q<16; q=q+2)
    {
      alpha = (vgrid[q+1].z - point.z)/(vgrid[q+1].z - vgrid[q].z);
      beta = 1 - alpha;      
      vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
      vgrid[q>>1] = vgrid[q];     
    }

  /* Lat Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(q = 0; q < 8; q=q+2)
    {
      vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
      vgrid[q>>1] = vgrid[q];
    }

  /* Phi Interpolation */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1.0 - alpha;
  for(q = 0; q < 4; q=q+2)
    {
      vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
      vgrid[q>>1] = vgrid[q];
    }

  /* Time Interpolation: */ 
  alpha = ((double) (time + 1)) - t;  
  beta = 1.0 - alpha; 
  vcomp[0] = alpha * vcomp[0] + beta * vcomp[1];

  /* Interpolated V*/
  *vint = vcomp[0];

  return 0;
}
