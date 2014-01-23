#ifndef SAMPLINGMC_H
#define SAMPLINGMC_H

#include "model.h"

using namespace std;

double new_beta (Model & model,
                 vector< vector<double> > & ensemble,
                 const size_t burn_in_step,
                 const size_t num_step,
                 const size_t max_num_step,
                 const size_t step_size,
                 const double a,
                 const double ,
                 const double C);
                 
double find_beta (Model & model, const vector<double> & Y, double , const double C, size_t);

size_t update_ensemble (Model & model,      // data and model
                        vector< vector<double> > & ensemble,
                        const double a);
                  
#endif
