/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *
 */
 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "data.h"
#include "int2str.h"
#include "model.h"
#include "rng.h"
#include "rosenbrock.h"

using namespace std;


// Constructor
Rosenbrock::Rosenbrock(): Model(2) {
	model_name = "rosenbrock";
	x1_upper_bound = 5;
	x1_lower_bound = -5;
	x2_upper_bound = 5;
	x2_lower_bound = -5;
	best_fit.resize(dim+1);
	best_fit[dim] = -numeric_limits<double>::infinity();
}

// Reparametrization (not used so far)
double Rosenbrock::reparametrize(const vector<double> & param, vector<double> & p) const
{
	for (size_t i = 0; i < dim; ++i) {
		p[i] = param[i];
	}
	return 1;
}

// Log Prior
// All the parameters have Uniform Prior
double Rosenbrock::LnPrior(const vector<double> & param) 
{
	vector<double> p(dim, 0);
	reparametrize(param, p);
	
	double lp = 0;
	
	if (p[0] < x1_lower_bound || p[0] > x1_upper_bound) {
		return -numeric_limits<double>::infinity();
	}
	if (p[1] < x2_lower_bound || p[1] > x2_upper_bound) {
		return -numeric_limits<double>::infinity();
	}
	
	lp = -log(x1_upper_bound - x1_lower_bound) - log(x2_upper_bound - x2_lower_bound);
	
	return lp;
}

// Log Likelihood
double Rosenbrock::LnLikelihood(const vector<double> & param)
{
	vector<double> p(dim, 0);
	reparametrize(param, p);
	
	double ll = - 5.0 * (p[1] - p[0]*p[0]) * (p[1] - p[0]*p[0]) - 0.05 * (1.0 - p[0]) * (1.0 - p[0]);
	
	return ll;
}
 
double Rosenbrock::LnDensity(const vector<double> & param) 
{
	if (beta.size() == 1) {
		return LnG(param);
	}
	else {
		return beta.back() * (LnPrior(param) + LnLikelihood(param)) + (1. - beta.back()) * LnG(param);
	}
}

double Rosenbrock::LnG(const vector<double> & parameter)
{
	vector<double> p(dim, 0);
	reparametrize(parameter, p);
	/*
	if (p[0] < x1_lower_bound || p[0] > x1_upper_bound) {
		return -numeric_limits<double>::infinity();
	}
	if (p[1] < x2_lower_bound || p[1] > x2_upper_bound) {
		return -numeric_limits<double>::infinity();
	}
	*/
	
	vector<double> temp(dim,0);
	
	
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			temp[i] += H[i][j] * (parameter[j]-m[j]);
		}
	}
	double exponent = 0;
	for (size_t i = 0; i < dim; ++i) {
		exponent += temp[i] * (parameter[i]-m[i]);
	}
	
	
	double logG = -0.5*(double)dim*log(2*M_PI) - 0.5*log_det_C - 0.5*exponent;
	
	return logG;
}

double Rosenbrock::LnImportance(const vector<double> & param)
{
	return LnPrior(param) + LnLikelihood(param) - LnG(param);
}

double Rosenbrock::LnPosterior(const vector<double> & param)
{
	return LnPrior(param) + LnLikelihood(param);
}

double Rosenbrock::fxx(const vector<double> & param) {
	vector<double> p(dim, 0);
	reparametrize(param, p);
	double lnlx = 20. * p[0] * (p[1] - p[0]*p[0]) + .1 * (1. - p[0]);
	double lnlxx = 20. * p[1] - 60. * p[0] - 0.1;
	double f = exp(LnPosterior(param));
	return f * lnlx * lnlx + f * lnlxx;
}

double Rosenbrock::fyy(const vector<double> & param) {
	vector<double> p(dim, 0);
	reparametrize(param, p);
	double lnly = -10. * (p[1] - p[0]*p[0]);
	double lnlyy = 10.;
	double f = exp(LnPosterior(param));
	return f * lnly * lnly + f * lnlyy;
}

// initialize the ensemble
void Rosenbrock::init(size_t ens_size, vector< vector<double> > & ensemble, double ini) const
{
	ensemble.resize(ens_size, vector<double>(dim, 0));
	for (size_t i = 0; i < ens_size; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			ensemble[i][j] = 0.01 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
		}
	}
}

Rosenbrock::~Rosenbrock() {
}


double rosenBruteIntegral(Rosenbrock & rosen, const size_t mesh_size_1, const size_t mesh_size_2) {
	double inc1 = (rosen.x1_upper_bound - rosen.x1_lower_bound) / static_cast<double>(mesh_size_1);        // increment in x1
	double inc2 = (rosen.x2_upper_bound - rosen.x2_lower_bound) / static_cast<double>(mesh_size_2);        // increment in x2
	double brute_integral = 0;
	double error = 0;
	vector<double> param(2);
	for (size_t i = 0; i < mesh_size_1; ++i) {
		for (size_t j = 0; j < mesh_size_2; ++j) {
			// use the center of the square, 2d rectangle rule?
			param[0] = rosen.x1_lower_bound + inc1 * (0.5 + i);
			param[1] = rosen.x2_lower_bound + inc2 * (0.5 + j);
			//cout << param[0] << "   " << param[1] << "   " << exp(rosen.LnPosterior(param)) << endl;
			error += 1./24. * rosen.fxx(param) * inc1 * inc1 + rosen.fyy(param) * inc2 * inc2;  // 2nd order term in Taylor series
			brute_integral += exp(rosen.LnPosterior(param));
		}
	}
	brute_integral *= inc1 * inc2;
	error          *= inc1 * inc2;
	//cout << "error = " << error << endl;
	return brute_integral;
}
