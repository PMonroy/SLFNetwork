#include "date.h"

int readparams(string nfileparameters);

/* PARAMATERS */

extern int domain_network;
extern int land_sea_norm; 
extern int velocity_field;
extern int verbose;       
extern double node_size;  
extern double particle_latspacing;
extern double start_depth;
extern double final_depth;
extern date start_date;
extern int tau; 
extern double int_step;
extern double v_sinking;

extern double degree_resolution;

extern double network_ll_lat;
extern double network_ll_lon;
extern double network_tr_lat;
extern double network_tr_lon;

extern double vfield_ll_lat;
extern double vfield_ll_lon;
extern double vfield_tr_lat;
extern double vfield_tr_lon;

/* CONFIGURATION PARAMETERS */
extern const double pig;
extern const double pig180;

/* Configuration parameters*/

extern const char output_dir[];
extern const char input_dir[];
// File names Outputs
extern const string postfixgrid;
extern const string postfixitracer;
extern const string postfixftracer;
extern const string postfixmatrix;

extern const char ncdir[];
extern date reference_date;
