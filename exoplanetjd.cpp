/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *
 */
 
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "arctan.h"
#include "covariance.h"
#include "data.h"
#include "exoplanet_hyperpara.h"
#include "exoplanet_init.h"
#include "exoplanetjd.h"
#include "factorial.h"
#include "int2str.h"
#include "keplers_eqn.h"
#include "linearalgebra.h"
#include "model.h"
#include "orbcrocor.h"
#include "sampling.h"

using namespace std;


size_t ExoplanetJD::count_distinction(void) {
	if (data.num_col == 3) {
		return 1;
	}
	if (data.num_col == 4) {
		vector<size_t> dist(data.num_row, 0);
		for (size_t i = 0; i < data.num_row; ++i) {
			dist[i] = static_cast<size_t>(data.data[i][3]); // data.data[][] is a double matrix.
		}
		
		vector<size_t> distinctions;
		distinctions.push_back(dist[0]);
		bool found;
		for (size_t i = 1; i < data.num_row; ++i) {
			found = 0;
			for (size_t j = 0; j < distinctions.size(); ++j) {
				if (dist[i] == distinctions[j]) {
					found = 1;
					break;
				}
			}
			if (found == 1) {
				continue;
			}
			else {
				distinctions.push_back(dist[i]);
			}
		}
		
		for (size_t i = 0; i < data.num_row; ++i) {
			found = 0;
			for (size_t j = 0; j < distinctions.size(); ++j) {
				if (dist[i] == distinctions[j]) {
					found = 1;
					data.data[i][3] = j;
					break;
				}
			}
			
			if (found != 1) {
				cerr << "Couldn't find a distinction index which should have been established!" << endl;
			}
		}
		return distinctions.size();
	}
	
	cerr << "To use ExoplanetJD, the data is at least required to have 3 or 4 columns." << endl;
	cerr << "The current data set Doesn't satisfy this least requirement." << endl;
	return 0;
}

// Constructor
ExoplanetJD::ExoplanetJD(Data & dt,                 // Data for likelihood calculation
                         size_t num_companion):     // Dimension of the model level
                         data(dt),
                         num_comp(num_companion){
	double begin_time;
	if (data.data_size != 0) {
		begin_time = data.data[0][0];
		for (size_t i = 0; i < data.data_size; ++i) {
			data.data[i][0] -= begin_time;
		}
	}
	model_name = "exop_" + data.data_name + "_mod_" + int2str(num_comp);
	num_d = count_distinction();
	if (num_d != 1) {
		model_name = model_name + "_d";
	}
	dim = 5 * num_comp + 2 * num_d;
	best_fit.resize(dim+1);
	best_fit[dim] = -1e300;
	LnPermutation = 0;
	for (size_t i = 1; i <= num_comp; ++i) {
		LnPermutation += log(static_cast<double>(i));
	}
	
	LnNorm = LnPermutation;
	bic_evi = 0.;
}

// Reparametrization
// input MCMC parameters
// output physical parameters
double ExoplanetJD::reparametrize(const vector<double> & param, vector<double> & p) const
{
	//param_0 = omega
	//param_1 = sqrt(amplitude) * sin(phi)
	//param_2 = sqrt(amplitude) * cos(phi)
	//param_3 = sqrt(eccentricity) * sin(varpi)
	//param_4 = sqrt(eccentricity) * cos(varpi)
	
	//This reparametrization has Jacobian 1.
	
	for (size_t i = 0; i < num_comp; ++i) {
		p[i*5+1] = param[i*5+0];
		p[i*5+0] = param[i*5+1] * param[i*5+1] + param[i*5+2] * param[i*5+2];
		p[i*5+3] = param[i*5+3] * param[i*5+3] + param[i*5+4] * param[i*5+4];
		p[i*5+2] = arctan(param[i*5+1], param[i*5+2]);  // The arctan takes two arguments 
		p[i*5+4] = arctan(param[i*5+3], param[i*5+4]);  // to have a range of [0, 2pi].
	}
	for (size_t i = 0; i < num_d; ++i) {
		p[num_comp*5 + i] = param[num_comp*5 + i];
	}
	for (size_t i = 0; i < num_d; ++i) {
		p[num_comp*5 + num_d + i] = param[num_comp*5 + num_d + i];
	}
	return 1; // return the jacobian
}

// Reparametrization
// input physical parameters
// output MCMC parameters
void ExoplanetJD::reverse_reparametrize(const vector<double> & p, vector<double> & param) const
{
	//param_0 = omega
	//param_1 = sqrt(amplitude) * sin(phi)
	//param_2 = sqrt(amplitude) * cos(phi)
	//param_3 = sqrt(eccentricity) * sin(varpi)
	//param_4 = sqrt(eccentricity) * cos(varpi)
	
	//This reparametrization has Jacobian 1.
	
	for (size_t i = 0; i < num_comp; ++i) {
		
		param[i*5+0] = p[i*5+1];
		param[i*5+1] = sqrt(p[i*5+0]) * sin(p[i*5+2]);
		param[i*5+2] = sqrt(p[i*5+0]) * cos(p[i*5+2]);
		param[i*5+3] = sqrt(p[i*5+3]) * sin(p[i*5+4]);
		param[i*5+4] = sqrt(p[i*5+3]) * cos(p[i*5+4]);
	}
	for (size_t i = 0; i < num_d; ++i) {
		param[num_comp*5 + i] = p[num_comp*5 + i] ;
	}
	for (size_t i = 0; i < num_d; ++i) {
		param[num_comp*5 + num_d + i] = p[num_comp*5 + num_d + i];
	}
}

// This routine checks if the orbits of two planets cross.
// If they do cross, return 1.
// If they do NOT cross, return 0.
bool ExoplanetJD::orbit_cross (const vector<double> & parameter) 
{  // parameters used in MCMC
	if (num_comp <= 1) {
		return 0;
	}

	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);

	double P_large, P_small;
	double a_large, a_small;
	double e_large, e_small;
	for (size_t i = 0; i < num_comp; ++i) {
		for (size_t j = i + 1; j < num_comp; ++j) {
			if (p[i*5+1] < p[j*5+1]) {
				P_large = 2 * M_PI / p[i*5+1];
				P_small = 2 * M_PI / p[j*5+1];
				a_large = pow(P_large*P_large, 1.0/3.0);
				a_small = pow(P_small*P_small, 1.0/3.0);
				e_large = p[i*5+3];
				e_small = p[j*5+3];
				if (a_small * (1 + e_small) > a_large * (1 - e_large)) {
					return 1;
				}
			}
			else {
				P_large = 2 * M_PI / p[j*5+1];
				P_small = 2 * M_PI / p[i*5+1];
				a_large = pow(P_large*P_large, 1.0/3.0);
				a_small = pow(P_small*P_small, 1.0/3.0);
				e_large = p[j*5+3];
				e_small = p[i*5+3];
				if (a_small * (1 + e_small) > a_large * (1 - e_large)) {
					return 1;
				}
			}
		}
	}
	return 0;
}

// If the periods of the companions are not in increasing order, return 1.
bool ExoplanetJD::exchange (const vector<double> & parameter) 
{
	if (num_comp <= 1) {
		return 0;
	}

	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);

	for (size_t i = 0; i < num_comp-1; ++i) {
		if (p[i*5+1] < p[(i+1)*5+1]) {
			return 1;
		}
	}
	return 0;
}

// If combinations of parameters are not physical, return 0.
bool ExoplanetJD::physical (const std::vector<double> & parameter)
{
	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	bool crossed = orbit_cross(parameter);
	//forbid orbit crossing
	if (crossed == 1) {
		return 0;  // 0 for unphysical
	}
	for (size_t i = 0; i < num_comp * 5; i=i+5) {
		if (p[i+3] < 0 || p[i+3] >=1) {
			return 0;
		}
		if (p[i+2] < 0) {
			return 0;
		}
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[5*num_comp+num_d+i] < 0) {
			return 0;
		}
	}
	return 1;
}

// the prior density function Uniform
double ExoplanetJD::LnPrior01 (const vector<double> & parameter)  // All Uniform Dummy Prior 
{

	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	bool crossed = orbit_cross(parameter);
	//forbid orbit crossing
	//if (crossed == 1) {
		//return -1e300;
	//}
	double lpr = 0;
	for (size_t i = 0; i < num_comp; ++i) {
		for (size_t j = 0; j < 5; ++j) {  // hard-wired FIVE because there are 5 orbital parameters.
			if (p[i*5+j] < hyper.uniOrbB[j][0] || p[i*5+j] > hyper.uniOrbB[j][1]) {
				return -numeric_limits<double>::infinity();
			}
			else {
				lpr += - log(hyper.uniOrbB[j][1] - hyper.uniOrbB[j][0]);
			}
		}
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+i] < hyper.uniV0B[0] || p[num_comp*5+i] > hyper.uniV0B[1]) {
			return -numeric_limits<double>::infinity();
		}
		else {
			lpr += - log(hyper.uniV0B[1] - hyper.uniV0B[0]);
		}
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+num_d+i] < hyper.uniJtqB[0] || p[num_comp*5+num_d+i] > hyper.uniJtqB[1]) {
			return -numeric_limits<double>::infinity();
		}
		else {
			lpr += - log(hyper.uniJtqB[1] - hyper.uniJtqB[0]);
		}
	}
	return lpr;
}	

// the prior density function More Complicated
double ExoplanetJD::LnPrior02 (const vector<double> & parameter) 
{
	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	//bool crossed = orbit_cross(parameter);
	if ( !physical(parameter) ) {  // unphysical
		return -numeric_limits<double>::infinity(); // small enough
	}
	if ( exchange(parameter) ) {
		return -numeric_limits<double>::infinity();
	}
	double lpr = 0;
	for (size_t i = 0; i < num_comp; ++i) {
		for (size_t j = 0; j < 5; ++j) {  // hard-wired FIVE because there are 5 orbital parameters.
			if (p[i*5+j] < hyper.uniOrbB[j][0] || p[i*5+j] > hyper.uniOrbB[j][1]) {
				return -numeric_limits<double>::infinity();
			}
		}
		lpr += - log(p[i*5+1] + hyper.JfOmega[0]);
		lpr += - log(p[i*5+0] + hyper.JfK[0]);
		lpr += (hyper.BetaEcc[0]-1)*log(p[i*5+3]) + (hyper.BetaEcc[1]-1)*log(1-p[i*5+3]);
		
		lpr += - log(log((hyper.uniOrbB[0][1]+hyper.JfK[0])/(hyper.uniOrbB[0][0]+hyper.JfK[0])));
		lpr += - log(log((hyper.uniOrbB[1][1]+hyper.JfOmega[0])/(hyper.uniOrbB[1][0]+hyper.JfOmega[0])));
		lpr += - log(hyper.uniOrbB[2][1] - hyper.uniOrbB[2][0]);
		lpr += log(factorial(static_cast<long>(hyper.BetaEcc[0]+hyper.BetaEcc[1]-1))) - log(factorial(static_cast<long>(hyper.BetaEcc[0]-1))) - log(factorial(static_cast<long>(hyper.BetaEcc[1]-1)));
		lpr += - log(hyper.uniOrbB[4][1] - hyper.uniOrbB[4][0]);
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+i] < hyper.uniV0B[0] || p[num_comp*5+i] > hyper.uniV0B[1]) {
			return -numeric_limits<double>::infinity();
		}
		lpr += - log(hyper.uniV0B[1] - hyper.uniV0B[0]);
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+num_d+i] < hyper.uniJtqB[0] || p[num_comp*5+num_d+i] > hyper.uniJtqB[1]) {
			return -numeric_limits<double>::infinity();
		}
		lpr += - log(p[num_comp*5+num_d+i] + hyper.JfJtq[0]);
		
		lpr += - log(log((hyper.uniJtqB[1]+hyper.JfJtq[0])/(hyper.uniJtqB[0]+hyper.JfJtq[0])));
	}
	lpr += LnNorm;
	return lpr;
}


// the logarithm of likelihood function
double ExoplanetJD::LnLikelihood(const vector<double> & parameter) 
{
	
	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	// If the eccentricity steps into absurd zones,
  // directly return a small value that means impossible.
	if ( !physical(parameter) ) {  // unphysical
		return -numeric_limits<double>::infinity(); // small enough
	}
	if ( exchange(parameter) ) {
		return -numeric_limits<double>::infinity();
	}

	double chi_sq      = 0.0;
	double ln_error_sq = 0.0;
	double rv_pred     = 0.0;
	double ln_jacobian = 0.0;
	double jitter_sq;
	
	for (size_t i = 0; i < data.data_size; ++i) {
		for (size_t j = 0; j < num_comp * 5; j=j+5) {
			rv_pred += rad_v_pred(data.data[i][0], p[j+0], p[j+1], p[j+2], p[j+3], p[j+4]);
			if (isnan(rv_pred) || isinf(rv_pred)) {
				cerr << "Keplers Equation Diverges." << endl;
				cerr << setprecision(17) << data.data[i][0] << "  " << p[j+0] << "  " << p[j+1] << "  " << p[j+2] << "  " << p[j+3] << "  " << p[j+4] << endl;
				cerr << "isnan = " << isnan(rv_pred) << ", isinf = " << isinf(rv_pred) << endl;
				return -numeric_limits<double>::infinity();
			}
		}
		if (data.num_col == 4) {
			rv_pred += p[num_comp*5+static_cast<size_t>(data.data[i][3])];
			jitter_sq = p[num_comp*5+num_d+static_cast<size_t>(data.data[i][3])];
		}
		if (data.num_col == 3) {
			rv_pred += p[num_comp*5];
			jitter_sq = p[num_comp*5+num_d];
		}
		if (data.num_col != 4 && data.num_col != 3) {
			cerr << "The data is required to have 3 or 4 columns." << endl;
			cerr << "Zero Likelihood will be returned." << endl;
			return -numeric_limits<double>::infinity();
		}
		chi_sq += (data.data[i][1] - rv_pred) * (data.data[i][1] - rv_pred) / (data.data[i][2] * data.data[i][2] + jitter_sq);  // the exponent
    
		ln_error_sq += log(2.0*M_PI*(data.data[i][2] * data.data[i][2] + jitter_sq)); // the normalization
		rv_pred = 0;
	}

	ln_jacobian = log(jacobian);
	
	double lnlikelihood = - 0.5 * ln_error_sq - 0.5 * chi_sq - ln_jacobian;
	if (lnlikelihood > best_fit.back()) {
		for (size_t k = 0; k < dim; ++k) {
			best_fit[k] = parameter[k];
		}
		best_fit[dim] = lnlikelihood;
	}
	return - 0.5 * ln_error_sq - 0.5 * chi_sq - ln_jacobian;
}

double ExoplanetJD::LnG(const vector<double> & parameter) {
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

double ExoplanetJD::LnDensity(const vector<double> & param)
{
	
	if (beta.size() == 1) {
		return LnG(param);
	}
	else {
		return beta.back() * (LnPrior02(param) + LnLikelihood(param)) + (1. - beta.back()) * LnG(param);
	}
}

double ExoplanetJD::LnImportance(const vector<double> & param)
{
	if ( !physical(param) ) {  // unphysical
		return -numeric_limits<double>::infinity();
	}
	if ( exchange(param) ) {
		return -numeric_limits<double>::infinity();
	}
	return LnPrior02(param) + LnLikelihood(param) - LnG(param) - bic_evi + static_cast<double>(num_comp)*log(4.0);
}

double ExoplanetJD::LnPosterior(const vector<double> & param)
{
	return LnPrior02(param) + LnLikelihood(param);
}

void ExoplanetJD::init(size_t ens_size, std::vector< std::vector<double> > & ensemble, double ini)
{
	bool fail = 0;
	srand(time(NULL));
	vector<double> param(dim);
	if (data.data_name.compare("data1") == 0 && num_d == 1) {
		ensemble.resize(ens_size, vector<double>(dim, 0));
		ExoplanetInit initialize(data.data_name, num_comp, hyper);
		reverse_reparametrize(initialize.init, param);
		for (size_t i = 0; i < ens_size; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				ensemble[i][j] = param[j] * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			}
		}
		
	}
	
}

void ExoplanetJD::init(size_t ens_size, std::vector< std::vector<double> > & ensemble, std::vector<double> & V, double ini)
{
	if (V.size() < dim) {
		cerr << "Not Enough Dimension in V" << endl;
	}
	ensemble.resize(ens_size, vector<double>(dim, 0));
	for (size_t i = 0; i < ens_size; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			ensemble[i][j] = V[j] * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
		}
	}
}

ExoplanetJD::~ExoplanetJD() {
}
