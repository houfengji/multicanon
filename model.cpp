/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  Feb 22, 2013
 *
 */
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "int2str.h"
#include "model.h"

using namespace std;

// Constructor:
// the dimension of the model's paramter space
Model::Model(size_t dimension):dim(dimension){

	time_label = int2str(static_cast<long>(time(NULL)));
	beta.push_back(0.);
	chain_evi.push_back(0.);
}

Model::Model(void) {
	time_label = int2str(static_cast<long>(time(NULL)));
	beta.push_back(0.);
	chain_evi.push_back(0.);
}


void Model::displayMeanCov(void) {
	cout << "Mean:" << endl << endl;
	for (size_t i = 0; i < m.size(); ++i) {
		cout << setprecision(16) << m[i] << endl;
	}
	cout << endl;
	cout << "Cov:" << endl << endl;
	for (size_t i = 0; i < C.size(); ++i) {
		for (size_t j = 0; j < C[0].size(); ++j) {
			cout << setprecision(16) << C[i][j] << "    ";
		}
		cout << endl;
	}
	cout << endl;
}

void Model::displayHessian(void) {
	cout << "Hessian:" << endl << endl;
	for (size_t i = 0; i < C.size(); ++i) {
		for (size_t j = 0; j < C[0].size(); ++j) {
			cout << setprecision(16) << H[i][j] << "    ";
		}
		cout << endl;
	}
	cout << endl;
}
