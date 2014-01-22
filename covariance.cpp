#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "acor.h"
#include "covariance.h"
#include "ensemblemean.h"
#include "model.h"
#include "sampling.h"

using namespace std;

int mean (const vector< vector<double> > & chain, vector<double> & m) {
	size_t L = chain.size();
	size_t d = chain[0].size();
	m.clear();
	m.resize(d, 0);
	
	for (size_t i = 0; i < d; ++i) {
		m[i] = 0;
		for (size_t k = 0; k < L; ++k) {
			m[i] += chain[k][i];
		}
		m[i] /= static_cast<double>(L);
	}
	
	return 1;
}

int covariance (const vector< vector<double> > & chain, vector< vector<double> > & cov) {
	
	size_t L = chain.size();
	size_t d = chain[0].size();
	cov.clear();
	cov.resize(d, vector<double>(d, 0.0));
	vector<double> m;
	mean(chain, m);
	for (size_t i = 0; i < d; ++i) {
		for (size_t j = 0; j < d; ++j) {
			cov[i][j] = 0.0;
			for (size_t k = 0; k < L; ++k) {
				cov[i][j] += (chain[k][i]-m[i]) * (chain[k][j]-m[j]);
			}
			cov[i][j] /= static_cast<double>(L);
		}
	}
	
	return 1;
}


int findMeanCov (ExoplanetJD & model,
                 vector< vector<double> > & ensemble,
                 const size_t per_cycle_steps,
                 const size_t num_burn,
                 const size_t num_cycle,
	               const double a) {
	
	int file_found = 0;
	string file_name_m, file_name_C;
	
	// If one already has the mean and covariance, one can load them from files.
	file_found = findMeanCovFromFile(model, file_name_m, file_name_C);
	if (file_found) {
		cout << "Mean and Hessian Files have been found!" << endl;
		return 1;
	}

	fstream chainout; // output the chain, optional bien sur.
	//chainout.open(("/export/bbq1/fh417/Gliese-581/chain_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out );

	for (size_t i = 0; i < num_burn; ++i) {
		for (size_t j = 0; j < per_cycle_steps; ++j) {
			sampling(model, ensemble, a);
		}
		if ( i % static_cast<long>(num_burn * 0.01) == 0 ) {
			cout << "burn-in " << i << " in " << num_burn << endl;
		}
	}
	
	
	
	vector< vector<double> > chain_ens;
	vector< vector<double> > chain_ens_mean;
	for (size_t i = 0; i < num_cycle; ++i) {
		for (size_t j = 0; j < per_cycle_steps; ++j) {
			sampling(model, ensemble, a);
		}
		if ( i % static_cast<long>(num_cycle * 0.01) == 0 ) {
			cout << "sample " << i << " in " << num_cycle << endl;
		}
		for (size_t k = 0; k < ensemble.size(); ++k) {
			chain_ens.push_back(ensemble[k]);
			for (size_t l = 0; l < model.dim; ++l) {
				chainout << "     " << setprecision(16) << ensemble[k][l];
			}
			chainout << endl;
		}
		chain_ens_mean.push_back( ensembleMean(ensemble) );
	}

	mean(chain_ens, model.m);
	covariance(chain_ens, model.C);
	
	/*
	vector<double> temp_chain;
	for (size_t i = 0; i < model.dim; ++i) {
		temp_chain.resize(0);
		for (size_t j = 0; j < chain_ens_mean.size(); ++j) {
			temp_chain.push_back(chain_ens_mean[j][i]);
		}
		double mean, sigma, tau;
		acor( mean, sigma, tau, temp_chain, temp_chain.size() );
		cout << "tau = " << tau << endl;
		model.C[i][i] *= tau;
	}
	*/

	
	fstream out;
	out.open(("h_cov_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	for (size_t i = 0; i < model.C.size(); ++i) {
		for (size_t j = 0; j < model.C[0].size(); ++j) {
			out << setprecision(17) << "    " << model.C[i][j];
		}
		out << endl;
	}
	out.close();
	out.open(("h_mean_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	for (size_t i = 0; i < model.m.size(); ++i) {
		out << setprecision(17) << "    " << model.m[i] << endl;
	}
	
	return 1;
}

int findMeanCov (ExoplanetJDPG & model,
                 vector< vector<double> > & ensemble,
                 const size_t per_cycle_steps,
                 const size_t num_burn,
                 const size_t num_cycle,
	               const double a) {
	
	int file_found = 0;
	string file_name_m, file_name_C;
	
	if ( model.model_name.compare( "exoppg_gliese_581pg_mod_3" ) == 0 ) {
		file_name_m = "./cov_mean/h_mean_exoppg_gliese_581pg_mod_3_1385495073.txt";
		file_name_C = "./cov_mean/h_cov_exoppg_gliese_581pg_mod_3_1385495073.txt";
	}
	if ( model.model_name.compare( "exoppg_gliese_581pg_mod_4" ) == 0 ) {
		file_name_m = "./cov_mean/h_mean_exoppg_gliese_581pg_mod_4_1385495520.txt";
		file_name_C = "./cov_mean/h_cov_exoppg_gliese_581pg_mod_4_1385495520.txt";
	}
	if ( model.model_name.compare( "exoppg_gliese_581pg_mod_5" ) == 0 ) {
		file_name_m = "./cov_mean/h_mean_exoppg_gliese_581pg_mod_5_1385501094.txt";
		file_name_C = "./cov_mean/h_cov_exoppg_gliese_581pg_mod_5_1385501094.txt";
	}
	if ( model.model_name.compare( "exoppg_gliese_581pg_mod_6" ) == 0 ) {
		file_name_m = "./cov_mean/h_mean_exoppg_gliese_581pg_mod_6_1385502576.txt";
		file_name_C = "./cov_mean/h_cov_exoppg_gliese_581pg_mod_6_1385502576.txt";
	}
	
	file_found = findMeanCovFromFile(model, file_name_m, file_name_C);
	if (file_found) {
		cout << "Mean and Hessian Files have been found!" << endl;
		return 1;
	}

	fstream chainout;
	//chainout.open(("/export/bbq1/fh417/Gliese-581/chain_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out );

	for (size_t i = 0; i < num_burn; ++i) {
		for (size_t j = 0; j < per_cycle_steps; ++j) {
			sampling(model, ensemble, a);
		}
		if ( i % static_cast<long>(num_burn * 0.01) == 0 ) {
			cout << "burn-in " << i << " in " << num_burn << endl;
		}
	}
	
	
	
	vector< vector<double> > chain_ens;
	vector< vector<double> > chain_ens_mean;
	for (size_t i = 0; i < num_cycle; ++i) {

		for (size_t j = 0; j < per_cycle_steps; ++j) {
			sampling(model, ensemble, a);
		}
		if ( i % static_cast<long>(num_cycle * 0.01) == 0 ) {
			cout << "sample " << i << " in " << num_cycle << endl;
		}
		for (size_t k = 0; k < ensemble.size(); ++k) {
			chain_ens.push_back(ensemble[k]);
			for (size_t l = 0; l < model.dim; ++l) {
				chainout << "     " << setprecision(16) << ensemble[k][l];
			}
			chainout << endl;
		}
		chain_ens_mean.push_back( ensembleMean(ensemble) );
	}

	mean(chain_ens, model.m);
	covariance(chain_ens, model.C);
	
	/*
	vector<double> temp_chain;
	for (size_t i = 0; i < model.dim; ++i) {
		temp_chain.resize(0);
		for (size_t j = 0; j < chain_ens_mean.size(); ++j) {
			temp_chain.push_back(chain_ens_mean[j][i]);
		}
		double mean, sigma, tau;
		acor( mean, sigma, tau, temp_chain, temp_chain.size() );
		cout << "tau = " << tau << endl;
		model.C[i][i] *= tau;
	}
	*/

	
	fstream out;
	out.open(("h_cov_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	for (size_t i = 0; i < model.C.size(); ++i) {
		for (size_t j = 0; j < model.C[0].size(); ++j) {
			out << setprecision(17) << "    " << model.C[i][j];
		}
		out << endl;
	}
	out.close();
	out.open(("h_mean_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	for (size_t i = 0; i < model.m.size(); ++i) {
		out << setprecision(17) << "    " << model.m[i] << endl;
	}
	
	return 1;
}

int findMeanCov (Model & model,
                 vector< vector<double> > & ensemble,
                 const size_t per_cycle_steps,
                 const size_t num_burn,
                 const size_t num_cycle,
	               const double a) {


	for (size_t i = 0; i < num_burn; ++i) {
		for (size_t j = 0; j < per_cycle_steps; ++j) {
			sampling(model, ensemble, a);
		}
		if ( i % static_cast<long>(num_burn * 0.01) == 0 ) {
			cout << "burn-in " << i << " in " << num_burn << endl;
		}
	}
	
	
	vector< vector<double> > chain_ens;
	vector< vector<double> > chain_ens_mean;
	for (size_t i = 0; i < num_cycle; ++i) {
		for (size_t j = 0; j < per_cycle_steps; ++j) {
			sampling(model, ensemble, a);
		}
		if ( i % static_cast<long>(num_cycle * 0.01) == 0 ) {
			cout << "sample " << i << " in " << num_cycle << endl;
		}
		for (size_t k = 0; k < ensemble.size(); ++k) {
			chain_ens.push_back(ensemble[k]);
		}
		chain_ens_mean.push_back( ensembleMean(ensemble) );
	}

	mean(chain_ens, model.m);
	covariance(chain_ens, model.C);
	
	/*
	vector<double> temp_chain;
	for (size_t i = 0; i < model.dim; ++i) {
		temp_chain.resize(0);
		for (size_t j = 0; j < chain_ens_mean.size(); ++j) {
			temp_chain.push_back(chain_ens_mean[j][i]);
		}
		double mean, sigma, tau;
		acor( mean, sigma, tau, temp_chain, temp_chain.size() );
		cout << "tau = " << tau << endl;
		model.C[i][i] *= tau;
	}
	*/

	
	fstream out;
	out.open(("h_cov_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	for (size_t i = 0; i < model.C.size(); ++i) {
		for (size_t j = 0; j < model.C[0].size(); ++j) {
			out << setprecision(17) << "    " << model.C[i][j];
		}
		out << endl;
	}
	out.close();
	out.open(("h_mean_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	for (size_t i = 0; i < model.m.size(); ++i) {
		out << setprecision(17) << "    " << model.m[i] << endl;
	}
	
	return 1;
}

int findMeanCovFromFile (Model & model, string & file_name_m, string & file_name_C) {
	size_t d = model.dim;
	model.m.clear();
	model.m.resize(d);
	model.C.clear();
	model.C.resize(d, vector<double>(d, 0.0));
	
	fstream in;
	in.open(file_name_m.c_str(), ios::in);
	if (in.fail()) {
		//cerr << "find m C: mean file: Failed to Open File " << file_name_m << endl;
		return 0;
	}
	for (size_t i = 0; i < d; ++i) {
		in >> model.m[i];
	}
	in.close();
	in.open(file_name_C.c_str(), ios::in);
	if (in.fail()) {
		//cerr << "find m C: cov file: Failed to Open File " << file_name_C << endl;
		return 0;
	}
	for (size_t i = 0; i < d; ++i) {
		for (size_t j = 0; j < d; ++j) {
			in >> model.C[i][j];
		}
	}
	in.close();
	return 1;
}

int find_m_file (Model & model, string & file_name_m) {
	size_t d = model.dim;
	model.m.clear();
	model.m.resize(d);
	model.C.clear();
	model.C.resize(d, vector<double>(d, 0.0));
	
	fstream in;
	in.open(file_name_m.c_str(), ios::in);
	if (in.fail()) {
		cerr << "find m C: mean file: Failed to Open File " << file_name_m << endl;
		return 0;
	}
	for (size_t i = 0; i < d; ++i) {
		in >> model.m[i];
	}
	in.close();
	return 1;
	
}

void enlarge_covariance (vector< vector<double> > & C, double df, double odf) {
	size_t dim = C.size();
	for (size_t i = 0; i < dim; ++i) {
		C[i][i] *= df;
	}
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			if (i != j) {
				C[i][j] *= odf;
			}
		}
	}
}
/*
int find_Hessian (Model & model) {
	
	if(model.Uplo != 'U' && model.Uplo != 'L') {
		cout << "wrong Uplo!" << endl;
		return 0;
	}
	
	if(model.Uplo == 'U') {
		vector< vector<double> > iU;
		vector< vector<double> > iUt;
		double A[model.invchC.size()];
		for (size_t i = 0; i < model.invchC.size(); ++i) {
			A[i] = model.invchC[i];
		}
		array_2_tri_matrix(&(model.Uplo), model.dim, A, iU);
		matrix_transpose(iU, iUt);
		matrix_multiplication(iU, iUt, model.H);
	}
	
	if(model.Uplo == 'L') {
		vector< vector<double> > iL;
		vector< vector<double> > iLt;
		double A[model.invchC.size()];
		for (size_t i = 0; i < model.invchC.size(); ++i) {
			A[i] = model.invchC[i];
		}
		array_2_tri_matrix(&(model.Uplo), model.dim, A, iL);
		matrix_transpose(iL, iLt);
		matrix_multiplication(iLt, iL, model.H);
	}
	
	fstream out;
	out.open("Hessian.txt" , ios::out);
	for (size_t i = 0; i < model.H.size(); ++i) {
		for (size_t j = 0; j < model.H.size(); ++j) {
			out << setprecision(15) << model.H[i][j] << "     ";
		}
		out << endl;
	}
	
	return 1;
}*/
