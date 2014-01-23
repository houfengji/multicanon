/*
 *  Fengji Hou
 *  New York University
 *  This cpp file contains routines that sample the mixture of
 *  posterior and gaussian, and finds beta increment (Delta beta_k).
 *  Nov 13, 2013
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "acor.h"
#include "ensemblemean.h"
#include "exception.h"
#include "int2str.h"
#include "mean.h"
#include "model.h"
#include "numzero.h"
#include "randn.h"
#include "rng.h"
#include "samplingmc.h"
#include "var.h"


#define TAUMAX  100              //  Compute tau directly only if tau < TAUMAX. Otherwise compute tau using the pairwise sum series.
#define WINMULT 10               //  Compute autocovariances up to lag s = WINMULT*TAU
#define MAXLAG  TAUMAX*WINMULT   //  The autocovariance array is double C[MAXLAG+1] so that C[s] makes sense for s = MAXLAG.
#define MINFAC   2               //  Stop and print an error message if the array is shorter than MINFAC * MAXLAG.

using namespace std;

// This routine finds beta_{k+1} given beta_k by sampling posterior^beta * gaussian^(1-beta).
// Return value is the acceptance ratio.
double new_beta (Model & model,
                 vector< vector<double> > & ensemble,
                 const size_t burn_in_step,
                 const size_t num_step,
                 const size_t max_num_step,
                 const size_t step_size,
                 const double a,        // tuning parameter of ensemble sampler
                 const double db_inc,   // increment of Delta beta
                 const double C) {      // the parameter that limits the increment of beta and bounds the error
	size_t ens_size = ensemble.size();
	size_t dim      = model.dim;
	const size_t max_subsample_step = num_step/(MINFAC*MAXLAG);
	
	
	if (model.dim != ensemble[0].size()) {
		cout << "WARNING: Model dim is not the same as Ensemble dim." << endl;
	}
	
	vector<double> Y;          // chain of log [ L(theta)pi(theta)/g(theta) ]
	vector<double> Emean;      // chain of the ensemble mean of Y
	vector<double> Emeansmall; // chain of the ensemble mean of some small number times Y
	vector<double> temp_ens;
	size_t accepted;
	
	size_t prac_step_size = step_size;
	double practical_a    = a;
	bool step_size_enlarged = 0;
	bool a_decreased = 0;
	
	vector<size_t> max_stuck_time(ens_size, 0);
	vector<size_t> current_stuck_time(ens_size, 0);
	vector< vector<double> > previous_ensemble(ens_size, vector<double>(dim, 0));
	size_t max_max_stuck_time = 0;
	bool changed = 0;
	for (size_t i = 0; i < burn_in_step; ++i) {
	
		for (size_t j = 0; j < ens_size; ++j) {
			for (size_t k = 0; k < dim; ++k) {
				previous_ensemble[j][k] = ensemble[j][k];
			}
		}
		
		accepted = 0;
		for (size_t j = 0; j < prac_step_size; ++j) {
			accepted += update_ensemble(model, ensemble, practical_a);
		}
		
		
		for (size_t j = 0; j < ens_size; ++j) {
			changed = 0;
			for (size_t k = 0; k < dim; ++k) {
				if (previous_ensemble[j][k] != ensemble[j][k]) {
					changed = 1;
					break;
				}
			}
			if ( !changed ) {
				++current_stuck_time[j];
				if (current_stuck_time[j] > max_stuck_time[j]) {
					max_stuck_time[j] = current_stuck_time[j];
				}
			}
			else {
				current_stuck_time[j] = 0;
			}
			/*
			if (current_stuck_time[j] > 100 && current_stuck_time[j]%100==0) {
				cout << "walker " << j << " has been stuck for " << current_stuck_time[j] << " times." << endl;
			}*/
		}
		
		max_max_stuck_time = *max_element(max_stuck_time.begin(), max_stuck_time.end());
		
		if ( i % (1+static_cast<long>(burn_in_step * 0.01)) == 0) {  // printing update on screen
			
			cout << "burn-in " << i << " in " << burn_in_step << " steps with step size " << prac_step_size << endl;
			cout << "maximum of the maximum of stuck time is " << max_max_stuck_time << endl;
			cout << "acceptance ratio is " << 1.0*accepted/ens_size/prac_step_size << endl;
		}
	}
	//cout << "Burn In Acc R: " << 1.0*accepted/ens_size/burn_in_step/step_size << endl;
	
	
	double logMean, logSigma, tau;

	accepted = 0;
	Emean.resize(0);
	Y.resize(0);
	for (size_t i = 0; i < max_num_step; ++i) {
		
		for (size_t j = 0; j < ens_size; ++j) {
			for (size_t k = 0; k < dim; ++k) {
				previous_ensemble[j][k] = ensemble[j][k];
			}
		}
		
		accepted = 0;
		for (size_t j = 0; j < prac_step_size; ++j) {
			accepted += update_ensemble(model, ensemble, practical_a);
		}
			
		for (size_t j = 0; j < ens_size; ++j) {
			changed = 0;
			for (size_t k = 0; k < dim; ++k) {
				if (previous_ensemble[j][k] != ensemble[j][k]) {
					changed = 1;
					break;
				}
			}
			if ( !changed ) {
				++current_stuck_time[j];
				if (current_stuck_time[j] > max_stuck_time[j]) {
					max_stuck_time[j] = current_stuck_time[j];
				}
			}
			else {
				current_stuck_time[j] = 0;
			}
		}
			
		max_max_stuck_time = *max_element(max_stuck_time.begin(), max_stuck_time.end());
			
		temp_ens.resize(0);
		for (size_t k = 0; k < ens_size; ++k) {
			temp_ens.push_back(model.LnImportance(ensemble[k]));
			Y.push_back(temp_ens[k]);
		}
		Emean.push_back(logEnsembleMeanLog(temp_ens));
		
		if ( i % (1+static_cast<long>(num_step * 0.01)) == 0 ) {  // printing update on screen
			cout << " sampling " << i << " in " << max_num_step << " steps with step size " << prac_step_size << endl;
			//cout << "maximum of the maximum of stuck time is " << max_max_stuck_time << endl;
			cout << "Emean size is " << Emean.size() << endl;
			cout << "acceptance ratio is " << 1.0*accepted/ens_size/prac_step_size << endl;
		}
		
		if ( i != 0 && Emean.size()%num_step == 0) {
			
			fstream out;
			//out.open(("./Emean/Emean_" + model.model_name+ "_" + model.time_label + "_bn_" + int2str((long)(model.beta.back()*10000),8) + "_" + int2str(Emean.size()) + ".txt").c_str(), ios::out);
			for (size_t j = 0; j < Emean.size(); ++j) {
				out << Emean[j] << endl;
			}
			out.close();
		
			double logMean, logSigma, tau;
			bool bad_tau = 0;
			try {
				logacorlog(logMean, logSigma, tau, Emean);
			}
			catch (Exception & e) {
				bad_tau = 1;
				fstream out;
				//out.open(("./Emean/Emean_" + model.model_name+ "_" + model.time_label + "_bb_" + int2str((long)(model.beta.back()*10000),8) + "_" + int2str(Emean.size()) + ".txt").c_str(), ios::out);
				for (size_t j = 0; j < Emean.size(); ++j) {
					out << Emean[j] << endl;
				}
				out.close();
			}
			
			vector<double> smallY(Y.size(),0);
			vector<double> smallEmean;
			for (size_t j = 0; j < Y.size(); ++j) {
				smallY[j] = Y[j] * db_inc;
			}
			smallEmean = logEnsembleMeanLog(smallY, ens_size);
			try {
				logacorlog(logMean, logSigma, tau, smallEmean);
			}
			catch (Exception & e) {
				bad_tau = 1;
				fstream out;
				//out.open(("./Emean/Emean_" + model.model_name+ "_" + model.time_label + "_bs_" + int2str((long)(model.beta.back()*10000),8) + "_" + int2str(Emean.size()) + ".txt").c_str(), ios::out);
				for (size_t j = 0; j < smallEmean.size(); ++j) {
					out << smallEmean[j] << endl;
				}
				out.close();
			}
			
			if ( !bad_tau ) {
				break;
			}
		}
	}

	//cout << "Sample Acc R: " << 1.0*accepted/ens_size/num_step/step_size << endl;
	
	//double tau_sub = static_cast<double>(max_max_stuck_time);
	
	
	
	double delta_beta = find_beta(model, Y, db_inc, C, ens_size);
	cout << "Delta Beta is " << delta_beta << endl;
	model.beta.push_back(model.beta.back()+delta_beta);
	
	bool bad_tau = 0;
	vector<double> Yb(Y.size());  // chain of log ( Importance Ratio ^ b )
	for (size_t i = 0; i < Y.size(); ++i) {
		Yb[i] = delta_beta * Y[i];
	}
	Emean = logEnsembleMeanLog(Yb, ens_size);
	try {
		logacorlog(logMean, logSigma, tau, Emean);
	}
	catch (Exception & e) {
		bad_tau = 1;
		cout << "Autocorrelation Error from new_beta apres obtenir Delta beta: " << endl;
		cout << e.ExceptionMessage() << endl;
		// output the chain for closer look
		fstream out;
		//out.open(("./Emean/Emean_" + model.model_name+ "_" + model.time_label + "_b_" + int2str(model.beta.back()*10000) + "_" + int2str(time(NULL)) + ".txt").c_str(), ios::out);
		for (size_t i = 0; i < Emean.size(); ++i) {
			out << Emean[i] << endl;
		}
		if (bad_tau) {
			tau = static_cast<double>(Emean.size()); // If error occurs when using acor, set tau to be the size of the chain.
		}
	}
	if (bad_tau) {
		model.tau.push_back(-tau);
	}
	else {
		model.tau.push_back(tau);
	}
	double logW = logMeanLog(Yb); // log ( mean ( exp(Yb) ) )
	double logWvar = logVarLog(Yb);
	double logR = 0.5 * ( logWvar + log(tau) - log(static_cast<double>(Yb.size())) ) - logW;
	model.chain_evi.push_back(model.chain_evi.back()+logW);
	model.chain_R.push_back(exp(logR));
	model.chain_C.push_back(C);
	
	return 1.0*accepted/ens_size/step_size;
}

double find_beta (Model & model, const vector<double> & Y, double db_inc, const double C, size_t ens_size) {

	const size_t max_subsample_step = Y.size()/(ens_size*MINFAC*MAXLAG);

	double db_min = 1.e-30;
	double db_max = 1. - model.beta.back();
	double db_new = 1.;
	
	vector<double> Yb(Y.size());
	vector<double> Emean;
	double logYbmean;
	double logYbvar;
	double logR;
	double logMean, logSigma, tau;
	
	// The scenario when there are 0's in the chain of Y
	size_t num_zero = numZeroLog(Y);
	cout << num_zero << " zeros in " << Y.size() << endl;
	if (num_zero > 0) {
		for (size_t i = 0; i < Y.size(); ++i) {
			Yb[i] = db_min * Y[i];
		}
		Emean = logEnsembleMeanLog(Yb, ens_size);
		try {
			logacorlog(logMean, logSigma, tau, Emean);
		}
		catch (Exception & e) {
			bool bad_tau = 1;
			cout << "Autocorrelation Error from find_beta zero: " << endl;
			cout << e.ExceptionMessage() << endl;
			if (bad_tau) {
				tau = static_cast<double>(Emean.size()); // If error occurs when using acor, set tau to be the size of the chain.
			}
		}
		double RLowLim = sqrt( static_cast<double>(num_zero) / static_cast<double>(Y.size() - num_zero) ) * sqrt(tau) / sqrt(static_cast<double>(Y.size()));
		if (RLowLim > C) {
			cout << "Delta beta_k = " << db_min << " beta_k = " << model.beta.back() << endl;
			model.chain_C.push_back(RLowLim);
			return db_min;
		}
	}
	
	// The scenario when 1 - beta_k is a good Delta beta_k
	for (size_t i = 0; i < Y.size(); ++i) {
		Yb[i] = db_max * Y[i];
	}
	Emean = logEnsembleMeanLog(Yb, ens_size);
	try {
		logacorlog(logMean, logSigma, tau, Emean);
	}
	catch (Exception & e) {
		cout << "Autocorrelation Error from find_beta 1-beta_k: " << endl;
		cout << e.ExceptionMessage() << endl;
		
		tau = static_cast<double>(Emean.size()); // If error occurs when using acor, set tau to be 100000 
	}
	logYbmean = logMeanLog(Yb);
	logYbvar = logVarLog(Yb);
	logR = 0.5 * ( logYbvar + log(tau) -log(static_cast<double>(Y.size())) ) - logYbmean;
	if (logR < log(C)) {
		cout << "Delta beta_k = " << db_max << " beta_k = " << model.beta.back() << endl;
		model.chain_C.push_back(C);
		return db_max;
	}
	
	for (size_t i = 0; i < Y.size(); ++i) {	
		Yb[i] = db_inc * Y[i];
	}
	Emean = logEnsembleMeanLog(Yb, ens_size);
	try {
		logacorlog(logMean, logSigma, tau, Emean);
	}
	catch (Exception & e) {
		cout << "Autocorrelation Error from new_beta db_inc: " << endl;
		cout << e.ExceptionMessage() << endl;
		if (db_max > 0.1) {
			db_max = 0.1;
		}
	}
	
	size_t count = 0;
	double good_tau = 0;
	double good_db = 0;
	bool tau_bad = 0;
	bool tau_once_good = 0;
	for (db_new = db_inc; db_new < db_max; db_new += db_inc) {
		for (size_t i = 0; i < Y.size(); ++i) {	
			Yb[i] = db_new * Y[i];
		}
		Emean = logEnsembleMeanLog(Yb, ens_size);
		try {
			logacorlog(logMean, logSigma, tau, Emean);
		}
		catch (Exception & e) {
			tau_bad = 1;
			bool bad_tau = 1;
			cout << "Autocorrelation Error from new_beta find_beta: " << endl;
			cout << e.ExceptionMessage() << endl;
			if (bad_tau) {
				tau = static_cast<double>(Emean.size()); // If error occurs when using acor, set tau to be the size of the chain.
			}
		}
		
		if (!tau_bad) {
			good_tau = tau;
			good_db  = db_new;
			tau_once_good = 1;
		}
		
		if (count != 0 && tau_once_good && tau_bad) {
			db_new = good_db;
			break;
		}
		logYbmean = logMeanLog(Yb);
		logYbvar = logVarLog(Yb);
		logR = 0.5 * ( logYbvar + log(tau) -log(static_cast<double>(Y.size())) ) - logYbmean;
		if (count == 0 && logR > log(C)) {
			db_inc = db_inc * 0.1;
			db_new = db_inc;
			count = 0;
			good_tau = 0;
			good_db = 0;
			tau_bad = 0;
			tau_once_good = 0;
			continue;
		}
		if (logR > log(C)) {
			db_new = db_new - db_inc;
			break;
		}
		
		if (count % 10 == 0) {
			cout << endl;
			cout << "Proposed Delta beta =  " << db_new << ", last beta = " << model.beta.back() << endl;
			cout << "log mean of Yb = " << logYbmean << endl;
			cout << "log var of Yb  = " << logYbvar << endl;
			cout << "tau            = " << tau << endl;
			cout << setprecision(10) << "R = " << exp(logR) << ", C = " << C << endl;
		}
		++count;
		
	}
	
	
	return db_new;
}
	
// the following routine updates the whole ensemble one step.
// Returned value is the number of moves accepted in the routine.
size_t update_ensemble (Model & model,      // data and model
                        vector< vector<double> > & ensemble,
                        const double a) {         // the tuning in the ensemble sampler

	size_t ens_size = ensemble.size();
	size_t dim      = model.dim;
	
	if (model.dim != ensemble[0].size()) {
		cout << "WARNING: Model dim is not the same as Ensemble dim." << endl;
	}
	
	vector<double> proposed_walker(dim, 0.0);
	size_t choose;
	double random, Z;
	double new_density, old_density;
	double accept;
	size_t accepted = 0;
	
	for (size_t k = 0; k < ens_size; ++k) {
		//choose a walker from the complementary ensemble which doesn't include walker_k
		int choose_fail = 0;
		do {
			choose = genrand_int32() % ens_size;
			
		} while (choose == k || choose == ens_size);
		
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
		
		new_density = model.LnDensity(proposed_walker);
		if (new_density < -1e200) {
			accept = 0;
		}
		else {
			old_density = model.LnDensity(ensemble[k]);
			if (new_density + (dim - 1.0) * log(Z) > old_density) {
				accept = 1;
			}
			else {
				accept = pow(Z, static_cast<int>(dim - 1.0)) * exp(new_density - old_density);
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
