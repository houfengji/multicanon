#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>
#include "exoplanetjd.h"
#include "factorial.h"
#include "int2str.h"
#include "orbcrocor.h"
#include "rng.h"

using namespace std;

double DensityOrbCroCor (ExoplanetJD & model,
                         const vector<double> & parameter) {
	vector<double> p(model.dim, 0);
	double jacobian = model.reparametrize(parameter, p);
	for (size_t i = 0; i < model.num_comp * 5; i=i+5) {
		if (p[i+3] < 0 || p[i+3] >=1) {
			return -numeric_limits<double>::infinity();;
		}
		if (p[i+2] < 0) {
			return -numeric_limits<double>::infinity();;
		}
	}
	for (size_t i = 0; i < model.num_d; ++i) {
		if (p[5*model.num_comp+model.num_d+i] < 0) {
			return -numeric_limits<double>::infinity();;
		}
	}
	if ( model.exchange(parameter) ) {
		return -numeric_limits<double>::infinity();
	}
	double lpr = 0;
	for (size_t i = 0; i < model.num_comp; ++i) {
		for (size_t j = 0; j < 5; ++j) {  // hard-wired FIVE because there are 5 orbital parameters.
			if (p[i*5+j] < model.hyper.uniOrbB[j][0] || p[i*5+j] > model.hyper.uniOrbB[j][1]) {
				return -numeric_limits<double>::infinity();
			}
		}
		lpr += - log(p[i*5+1] + model.hyper.JfOmega[0]);
		lpr += - log(p[i*5+0] + model.hyper.JfK[0]);
		lpr += (model.hyper.BetaEcc[0]-1)*log(p[i*5+3]) + (model.hyper.BetaEcc[1]-1)*log(1-p[i*5+3]);
		
		lpr += - log(log((model.hyper.uniOrbB[0][1]+model.hyper.JfK[0])/(model.hyper.uniOrbB[0][0]+model.hyper.JfK[0])));
		lpr += - log(log((model.hyper.uniOrbB[1][1]+model.hyper.JfOmega[0])/(model.hyper.uniOrbB[1][0]+model.hyper.JfOmega[0])));
		lpr += - log(model.hyper.uniOrbB[2][1] - model.hyper.uniOrbB[2][0]);
		lpr += log(factorial(static_cast<long>(model.hyper.BetaEcc[0]+model.hyper.BetaEcc[1]-1))) - log(factorial(static_cast<long>(model.hyper.BetaEcc[0]-1))) - log(factorial(static_cast<long>(model.hyper.BetaEcc[1]-1)));
		lpr += - log(model.hyper.uniOrbB[4][1] - model.hyper.uniOrbB[4][0]);
	}
	for (size_t i = 0; i < model.num_d; ++i) {
		if (p[model.num_comp*5+i] < model.hyper.uniV0B[0] || p[model.num_comp*5+i] > model.hyper.uniV0B[1]) {
			return -numeric_limits<double>::infinity();
		}
		lpr += - log(model.hyper.uniV0B[1] - model.hyper.uniV0B[0]);
	}
	for (size_t i = 0; i < model.num_d; ++i) {
		if (p[model.num_comp*5+model.num_d+i] < model.hyper.uniJtqB[0] || p[model.num_comp*5+model.num_d+i] > model.hyper.uniJtqB[1]) {
			return -numeric_limits<double>::infinity();
		}
		lpr += - log(p[model.num_comp*5+model.num_d+i] + model.hyper.JfJtq[0]);
		
		lpr += - log(log((model.hyper.uniJtqB[1]+model.hyper.JfJtq[0])/(model.hyper.uniJtqB[0]+model.hyper.JfJtq[0])));
	}
	return lpr;
}
                         

size_t SamplingOrbCroCor (ExoplanetJD & model,
                          vector< vector<double> > & ensemble,
                          const double a) {         // the tuning in the ensemble sampler

	size_t ens_size = ensemble.size();
	size_t dim      = ensemble[0].size();
	
	vector<double> proposed_walker(dim, 0.0);
	size_t choose;
	double random, Z;
	double new_density, old_density;
	double accept;
	size_t accepted = 0;
	
	for (size_t k = 0; k < ens_size; ++k) {
		//choose a walker from the complementary ensemble which doesn't include walker_k
		do {
			choose = genrand_int32() % ens_size;
		} while (choose == k);
		
		random = genrand_real2();
		//Z is drawn from a distribution satisfying g(z)=g(1/z)/z.
		//The distribution recommanded in Goodman and Weare's paper is used here.
		//To sample this distribution, direct sampling is the easiest.
		Z = ((a - 1.0) * random + 1.0) * ((a - 1.0) * random + 1.0) / a;
		
		//proposal based on stretch move
		for (size_t j = 0; j < dim; ++j) {
			//X_j(t+1) = Y_j(t) + Z * (X_j(t) - Y_j(t))
			//where Y belongs to the complementary ensemble
			proposed_walker[j] = ensemble[choose][j] + Z * (ensemble[k][j] - ensemble[choose][j]);
			
		}
		
		new_density = DensityOrbCroCor(model, proposed_walker);
		if (new_density == -numeric_limits<double>::infinity()) {
			accept = 0;
		}
		else {
			old_density = DensityOrbCroCor(model, ensemble[k]);
			if (new_density + (dim - 1.0) * log(Z) >= old_density) {
				accept = 1;
			}
			else {
				accept = pow(Z, static_cast<double>(dim - 1.0)) * exp(new_density - old_density);
			}
		}

		//accept or reject based on accept
		random = genrand_real2();
		if (accept > random) {
			for (size_t j = 0; j < dim; ++j) {
				ensemble[k][j] = proposed_walker[j];
			}
			//ensemble[k] = proposed_walker;
			accepted += 1;
		}
	}
	return accepted;
}

int OrbCroCor (ExoplanetJD & model, 
               vector< vector<double> > & ensemble,
               const size_t per_cycle_steps,
               const size_t num_burn,
               const size_t num_cycle,
	             const double a) {

	if (model.num_comp == 1) {
		cout << "Num of Comp is 1." << endl;
		return 1;
	}
	
	bool norm_found = 0;
	
	if (norm_found) {
		cout << "Normalization of Prior found!" << endl;
		return 1;
	}
	
	for (size_t i = 0; i < num_burn; ++i) {
		for (size_t j = 0; j < per_cycle_steps; ++j) {
			SamplingOrbCroCor(model, ensemble, a);
		}
		if ((double)rand()/(double)RAND_MAX < 0.01) {
			cout << "burn-in " << i << " in " << num_burn << endl;
		}
	}
	
	
	double num_uncross = 0;
	double total = 0;
	for (size_t i = 0; i < num_cycle; ++i) {
		for (size_t j = 0; j < per_cycle_steps; ++j) {
			SamplingOrbCroCor(model, ensemble, a);
		}
		if ((double)rand()/(double)RAND_MAX < 0.01) {
			cout << "orbital crossing correstion, sample " << i << " in " << num_cycle << endl;
		}
		for (size_t k = 0; k < ensemble.size(); ++k) {
			if (model.orbit_cross(ensemble[k]) == 1) {
				total = total + 1.0;
			}
			else {
				num_uncross = num_uncross + 1.0;
				total = total + 1.0;
			}
		}
	}
	
	model.LnNorm += log(total/num_uncross);
	fstream out;
	out.open(("Cross_discount_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	out << "num burn = " << num_burn << endl;
	out << "num cycle = " << num_cycle << endl;
	out << "step size = " << per_cycle_steps << endl;
	out << "a = " << a << endl;
	out << "ensemble size = " << ensemble.size() << endl;
	out << "K_min = " << model.hyper.uniOrbB[0][0] << " m/s" << endl;
	out << "K_max = " << model.hyper.uniOrbB[0][1] << " m/s" << endl;
	out << "omega_min = " << model.hyper.uniOrbB[1][0] << " rad/d" << endl;
	out << "omega_max = " << model.hyper.uniOrbB[1][1] << " rad/d" << endl;
	out << setprecision(15) << num_uncross << " un-crossed in " << total << endl;
	out << "discount = " << setprecision(15) << log(total/num_uncross) << endl;
	
	return 1;
}
