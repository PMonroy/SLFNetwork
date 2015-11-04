#include "date.h"

enum  enum_parameters { /* 0 */ NETWORKDOMAIN,
			/* 1 */ VFLOW,
			/* 2 */ VERBOSE,
			/* 3 */ NODESIZE,
			/* 4 */ PARTICLESPACING,
			/* 5 */ STARTDEPTH,
			/* 6 */ FINALDEPTH,
			/* 7 */ STARTDATE,
			/* 8 */ INTSTEP,
			/* 9  */ VSINK,
			NPARAMETERS
};

int readparams(string nfileparameters);
