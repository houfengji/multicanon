/*
 *  main.cpp
 */
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "acor.h"
#include "covariance.h"
#include "data.h"
#include "exception.h"
#include "exoplanetjd.h"
#include "int2str.h"
#include "linearalgebra.h"
#include "mean.h"
#include "model.h"
#include "orbcrocor.h"
#include "record.h"
#include "rng.h"
#include "rosenbrock.h"
#include "sampling.h"
#include "samplingmc.h"
#include "var.h"


using namespace std;

int read_beta (Model & model, const string & filename, const int head);

int main(void) {
	time_t begin = time(NULL);
	init_genrand(begin);
	time_t end;
	
	
	string data_file_name("data1"); size_t num_comp = 1; size_t ens_size = 100; size_t step_size = 5;
	Data data(data_file_name);
	ExoplanetJD model(data, num_comp);
	cout << model.model_name << endl;
	vector< vector<double> > ensemble;
	double ini = 0.00001;
	model.init(ens_size, ensemble, ini);
	cout << "Normalizatin = " << exp(model.LnNorm) << endl;
	
	step_size = 10;
	size_t num_burn = (1+num_comp)*4000;
	size_t num_cycle = (1+num_comp)*20000;
	double aa = 2;
	//cout << " Sample LnPosterior = " << model.LnPosterior(ensemble[0]) << endl;
	model.init(ens_size, ensemble, ini);
	OrbCroCor(model, ensemble, step_size, num_burn, num_cycle, aa);
	model.init(ens_size, ensemble, ini);
	findMeanCov(model, ensemble, step_size, num_burn, num_cycle, aa);
	//cout << model.LnNorm - model.LnPermutation << endl;
	try {
		matrix_inverse_cholesky(model.C, model.H);
	}
	catch (Exception & e) {
		cout << "Inverse of Covariance Failed:" << endl;
		cout << e.ExceptionMessage() << endl;
	}
	//model.displayHessian();
	
	try {
		model.log_det_C = matrix_log_determinant(model.C);
	}
	catch (Exception & e) {
		cout << "Determinant of Covariance Calculation Failed:" << endl;
		cout << e.ExceptionMessage() << endl;
	}
	vector<double> maximum_likelihood(model.dim, 0.);
	for (size_t i = 0; i < model.dim; ++i) {
		maximum_likelihood[i] = model.best_fit[i];
	}
	model.bic_evi =  0.5 * model.dim * log(2.0*M_PI) + 0.5 * model.log_det_C + model.LnPosterior(maximum_likelihood);
	
	
	size_t burn_in_step = 2000*(num_comp+1);
	size_t num_step = 10000*(num_comp+1);
	step_size = 10;
	double a = 1.5;                             // tuning parameter of sampler
	double C = 0.01;                            // beta increase variance control
	double C_change = 0.6;
	double C_later = 0.01;                       // beta increase variance control
	double a_later = 1.2;
	double db_inc = 1.e-3;                      // Delta beta accuracy
	
	ini = 0.001;
	
	fstream out;
	out.open(("results_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::app | ios::out);
	out << "burn-in steps: " << burn_in_step << endl;
	out << "num of steps:  " << num_step << endl;
	out << "ensemble size: " << ens_size << endl;
	out << "step size:     " << step_size << endl;
	out << "C:             " << setprecision(16) << C << endl;
	out << "a:             " << a << endl;
	out << "cutoff:        " << C_change << endl;
	out << "C_later:       " << C_later << endl;
	out << "a_labter:      " << a_later << endl;
	out << "K_min:         " << model.hyper.uniOrbB[0][0] << endl;
	out << "K_max:         " << model.hyper.uniOrbB[0][1] << endl;
	out << "omega_min:     " << model.hyper.uniOrbB[1][0] << endl;
	out << "omega_max:     " << model.hyper.uniOrbB[1][1] << endl;
	out << "BIC evidence:  " << model.bic_evi << endl;
	
	size_t i = 0;
	double accptr = 0;
	while(  model.beta.back() < 1)  {
		++i;
		
		if (model.beta.size() == 1) {
			model.init(ens_size, ensemble, model.m, ini); // initializing ensemble
			// longer burn-in for better mixing
			new_beta(model, ensemble, burn_in_step*5, num_step, num_step*10, step_size, a, db_inc, C); 
			model.init(ens_size, ensemble, ini); // initializing ensemble
		}
		else {
			if (model.beta.size() == 2) {
				// Due to the re-initialization, longer burn-in is needed for the 2nd Delta beta.
				new_beta(model, ensemble, burn_in_step*5/2, num_step, num_step*10, step_size, a, db_inc, C); 
			}
			else {
				if (model.beta.back() < C_change) {
					accptr = new_beta(model, ensemble, burn_in_step, num_step, num_step*10, step_size, a, db_inc, C); 
				}
				else {
					accptr = new_beta(model, ensemble, burn_in_step*10, num_step, num_step*10, step_size, a_later, db_inc, C_later);
				}
			}
		}
		record(model, out);
	}
	out << "Number of Inc : " << model.chain_evi.size()-1 << endl;
	
	out << "Log Evidence  : " << model.bic_evi + model.chain_evi.back() << endl;
	out << "Relative Error: " << sqrt(model.chain_evi.size()-1)*C << endl;
	
	end = time(NULL);
	out << "Time Used:      " << end-begin << " seconds!" << endl;
	
	
	
	//fstream aout;
	//aout.close();
	//aout.open(("resultmc_"+model.model_name+"_"+"2013"+".txt").c_str(), ios::app | ios::out);
	//aout << setprecision(16) << model.chain_evi.back() << "     " << logSTD << "   " << model.chain_evi.size()-1 << "    " << C << "     " << ens_size*num_step << endl;
	//aout.close();
	
}

int read_beta (Model & model, const string & filename, const int head) {
	fstream in;
	in.open(filename.c_str(), ios::in);
	if (in.fail()) {
		return 0;
	}
	string temporary; // temporary string used to store the line read by getline()
	for (int i = 0; i < head; ++i) {
		getline( in, temporary );
	}
	double temp[6];
	do {
		if (in >> temp[0] >> temp[1] >> temp[2] >> temp[3] >> temp[4] >> temp[5] ) {
			model.beta.push_back(temp[0]);
			model.chain_evi.push_back(temp[1]);
			model.chain_C.push_back(temp[3]);
			model.chain_R.push_back(temp[4]);
			model.tau.push_back(temp[5]);
		}
	}while(!in.eof());
	return 1;
}
