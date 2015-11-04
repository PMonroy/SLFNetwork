#ifndef VELOCITY
#define VELOCITY

#include "date.h"
#include "vectorXYZ.h"

/* Transformations: */
#define PI 3.1415926535898
#define RADS(degree) (PI/180.0)*(degree)
#define DEGREE(rads) (180.0/PI)*(rads)
#define MU(x)  log(fabs((1.0/cos(x))+tan(x))) // x must be in radians. Result value in radians
#define THETA(x) 1.56864 -1.99805 * atan(exp(-1.0*(x))) // x must be in radians. Result value in radians

/* Structures: Esta estructura tiene que ser interna a velocity.cpp */
struct  vectorIJK {
    int i;
    int j;
    int k[2][2];
};

/* FUNCTIONS */
int LoadVelocityGrid(date reference_date, char velocitydir[]);
void FreeMemoryVelocityGrid();
int LoadVelocities(date start_date, int tau, char velocitydir[]);// This function read the 2D velocity field since start_date to start_date+tau in a constant depth (layer_index)
void FreeMemoryVelocities(int tau);

int GetLonIndex(double latitude);
int GetLatIndex(double longitude);
int GetIndices(unsigned long time, vectorXYZ point, vectorIJK *index);
int GetVelocity(double t,vectorXYZ point, vectorXYZ *vint);

/* Global Variables*/
extern int nlon, nlat, ndepth, ntime;
extern double *vgrid_lon, *vgrid_lat, ****vgrid_depth;;
extern vectorXYZ ****vfield;
extern int **land_mask;
extern double **bathymetry;


#define SECONDS_DAY 86400.0

#endif
