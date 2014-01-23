#ifndef ORBCROCOR_H
#define ORBCROCOR_H

#include <vector>
#include "exoplanetjd.h"

using namespace std;

int OrbCroCor (ExoplanetJD & model, 
               vector< vector<double> > & ensemble,
               const size_t per_cycle_steps,
               const size_t num_burn,
               const size_t num_cycle,
	             const double a);
	                             
#endif
