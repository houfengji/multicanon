#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include "mean.h"

using namespace std;

vector<double> ensembleMean(const vector< vector<double> > & ensemble) {
	size_t ens_size = ensemble.size();
	size_t dim   = ensemble[0].size();
	vector<double> Emean(dim, 0);
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = 0; j < ens_size; ++j) {
			Emean[i] += ensemble[j][i];
		}
		Emean[i] /= static_cast<double>(ens_size);
	}
	return Emean;
}

vector<double> ensembleMean(const vector<double> & chain, const size_t ens_size) {
	
	size_t length = chain.size() / ens_size;
	vector<double> mean_chain;
	if (length * ens_size != chain.size()) {
		cerr << "Chain length cannot be completely divided by ensemble size!" << endl;
		return mean_chain;
	}
	vector<double> ensemble;
	for (size_t i = 0; i < length; ++i) {
		ensemble.resize(0);
		for (size_t k = 0; k < ens_size; ++k) {
			ensemble.push_back(chain[i*ens_size+k]);
		}
		mean_chain.push_back( mean(ensemble) );
	}
	return mean_chain;
}

double logEnsembleMeanLog(const vector<double> & ensemble) {
	size_t ens_size = ensemble.size();
	double max_ens = *max_element(ensemble.begin(), ensemble.end());
	if (max_ens == -numeric_limits<double>::infinity()) {
		return -numeric_limits<double>::infinity();
	}
	double Emean = 0;
	for (size_t i = 0; i < ens_size; ++i) {
		Emean += exp(ensemble[i] - max_ens);
	}
	Emean /= static_cast<double>(ens_size);
	return log(Emean) + max_ens;
}

vector<double> logEnsembleMeanLog(const vector<double> & chain, const size_t ens_size) {
	
	size_t length = chain.size() / ens_size;
	vector<double> mean_chain;
	if (length * ens_size != chain.size()) {
		cerr << "Chain length cannot be completely divided by ensemble size!" << endl;
		return mean_chain;
	}
	vector<double> ensemble;
	for (size_t i = 0; i < length; ++i) {
		ensemble.resize(0);
		for (size_t k = 0; k < ens_size; ++k) {
			ensemble.push_back(chain[i*ens_size+k]);
		}
		mean_chain.push_back( logMeanLog(ensemble) );
	}
	return mean_chain;
}
