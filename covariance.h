#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <vector>
#include "model.h"
#include "exoplanetjd.h"
using namespace std;

int mean (const vector< vector<double> > & chain, vector<double> & m);
int covariance (const vector< vector<double> > & chain, vector< vector<double> > & cov);

int findMeanCov (ExoplanetJD & model,
                 vector< vector<double> > & ensemble, 
                 const size_t per_cycle_steps,
                 const size_t num_burn,
                 const size_t num_cycle,
	               const double a);
	               
int findMeanCov (Model & model,
                 vector< vector<double> > & ensemble, 
                 const size_t per_cycle_steps,
                 const size_t num_burn,
                 const size_t num_cycle,
	               const double a);
	            
int findMeanCovFromFile (Model & model, string & m, string & C);
int find_m_file (Model & model, string & file_name_m);
void enlarge_covariance (vector< vector<double> > & C, double df, double odf) ;

#endif
