#include "vectorXYZ.h"
/* To use file I/O */
#include <iostream>
#include <fstream>
/* Standard libs */
#include <cstdlib>
#include <cstdio>
/* I dont know */
#include <cmath>

using namespace std;

void MakeVTKStructuredGrid_Vfield(double *lon, double *lat, double ***depth, vectorXYZ ***vfield, int ni, int nj, int nk, string filename)
{
  int i,j,k;
  ofstream ofile(filename.c_str());

  ofile<<"# vtk DataFile Version 3.0"<<endl;
  ofile<<"Complete vector field of ROMS Benguela"<<endl; 
  ofile<<"ASCII"<<endl;
  ofile<<"DATASET STRUCTURED_GRID"<<endl;
  ofile<<"DIMENSIONS "<<ni<<" "<<nj<<" "<<nk<<endl;
  ofile<<"POINTS "<<ni*nj*nk<<" float"<<endl;
  for(k=0;k<nk;k++) 
    {
      for(j=0;j<nj;j++) 
	{
	  for(i=0;i<ni;i++)
	    {		
	      ofile << lon[i] <<" "<< lat[j]<<" "<< depth[i][j][k]/1000.0<<endl;
	    }
	}
    }
  ofile<<endl;
  ofile<<"POINT_DATA "<<ni*nj*nk<<endl;
  ofile<<"SCALARS depth float"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(k=0;k<nk;k++) 
    {
      for(j=0;j<nj;j++) 
	{
	  for(i=0;i<ni;i++)
	    {		
	      ofile<<depth[i][j][k]<<endl;
	    }
	}
    }

  ofile<<"VECTORS velocity float"<<endl;
  for(k=0;k<nk;k++) 
    {
      for(j=0;j<nj;j++) 
	{
	  for(i=0;i<ni;i++)
	    {		
	      ofile<<vfield[i][j][k].x<<" "<<vfield[i][j][k].y<<" "<<vfield[i][j][k].z<<endl;
	    }
	}
    }

  ofile.close();
}
void MakeVTKStructuredGrid2D(double *lon, double *lat, double depth, double **A, int ni, int nj, string filename)
{
  int i,j,k;
  ofstream ofile(filename.c_str());

  ofile<<"# vtk DataFile Version 3.0"<<endl;
  ofile<<"Complete vector field of ROMS Benguela"<<endl; 
  ofile<<"ASCII"<<endl;
  ofile<<"DATASET STRUCTURED_GRID"<<endl;
  ofile<<"DIMENSIONS "<<ni<<" "<<nj<<" "<<1<<endl;
  ofile<<"POINTS "<<ni*nj<<" float"<<endl;
  for(j=0;j<nj;j++) 
    {
      for(i=0;i<ni;i++)
	{		
	  ofile << lon[i] <<" "<< lat[j]<<" "<< A[i][j]/1000.0<<endl;
	}
    }
  ofile<<endl;
  ofile<<"POINT_DATA "<<ni*nj<<endl;
  ofile<<"SCALARS bathymetry float"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(j=0;j<nj;j++) 
    {
      for(i=0;i<ni;i++)
	{		
	  ofile<<A[i][j]<<endl;
	}
    }

  ofile.close();
}
