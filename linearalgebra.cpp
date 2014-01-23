/*
 *   Fengji Hou
 *   fh417@nyu.edu
 *   New York University
 *   This cpp file contains the LA related functions.
 *
 */

//#include "/home/fh417/Projects_Local/CLAPACK-3.2.1/f2c.h"
//#include "/home/fh417/Projects_Local/CLAPACK-3.2.1/clapack.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <vector>
#include "exception.h"
#include "linearalgebra.h"

using namespace std;

	
bool symmetric (const vector< vector<double> > & C) {
	if (C.size() != C[0].size()) {
		cerr << "num of rows != num of columns" << endl;
		return 0;
	}
	if (C.size() == 0) {
		cerr << "nothing in the matrix" << endl;
		return 0;
	}
	if (C.size() == 1) {
		cout << "scalar matrix" << endl;
		return 1;
	}
	for (size_t i = 0; i < C.size(); ++i) {
		for (size_t j = i+1; j < C.size(); ++j) {
			if (fabs(C[i][j] - C[j][i]) > fabs(C[i][j])*1.e-17) {
				cerr << "not symmetric" << endl;
				return 0;
			}
		}
	}
	return 1;
}

// dumb version of matrix multiplication
int matrix_multiplication ( vector< vector<double> > & A,
                            vector< vector<double> > & B,
                            vector< vector<double> > & C) {

	if ( A.size() == 0 || B.size() == 0) {
		cerr << "matrix_multiplication: empty matrix." << endl;
		return 0;
	}
	if ( A[0].size() != B.size() ) {
		cerr << "matrix_multiplication: num col of A != num row of B!" << endl;
		return -1;
	}
	
	C.clear();
	C.resize( A.size(), vector<double>(B[0].size(), 0) );
	
	for (size_t i = 0; i < C.size(); ++i) {
		for (size_t j = 0; j < C[0].size(); ++j) {
			for (size_t k = 0; k < A[0].size(); ++k) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	
	return 1;
}

int matrix_transpose( vector< vector<double> > & A,
                      vector< vector<double> > & At) {

	if ( A.size() == 0) {
		cerr << "matrix_transpose: empty matrix." << endl;
		return 0;
	}
	
	At.clear();
	At.resize( A[0].size(), vector<double>(A.size(), 0.) );
	
	for (size_t i = 0; i < At.size(); ++i) {
		for (size_t j = 0; j < At[0].size(); ++j) {
			At[i][j] = A[j][i];
		}
	}
	return 1;
}

int matrix_cholesky( const vector< vector<double> > & A,
                     vector< vector<double> > & R,
                     char uplo) {

	if (A.size() != A[0].size()) {
		cerr << "matrix_cholesky: Matrix Not Square!" << endl;
		return 0;
	}
	
	long N = A.size();    // dimension of the matrix
	R.clear();
	R.resize( N, vector<double>(N, 0.) );
	
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = i; j < N; ++j) {
			if ( i != j && A[i][j] != A[j][i] ) {
				cerr << "matrix_cholesky: Matrix Not Symmetric!" << endl;
				return 0;
			}
			R[i][j] = A[i][j];
		}
	}
	
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = i+1; j < N; ++j) {
			for (size_t k = j; k < N; ++k) {
				R[j][k] = R[j][k] - R[i][k] * R[i][j] / R[i][i];
			}
		}
		double Rii = R[i][i];
		if (Rii < 0.) {
			cerr << "matrix_cholesky: Matrix Not Positive Definite!" << endl;
			return 0;
		}
		for (size_t j = i; j < N; ++j) {	
			R[i][j] = R[i][j] / sqrt(Rii);
		}
	}
	
	if (uplo == 'U' || uplo == 'u') {
		return 1;
	}
	if (uplo == 'L' || uplo == 'l') {
		for (size_t i = 0; i < N; ++i) {
			for (size_t j = i+1; j < N; ++j) {
				R[j][i] = R[i][j];
				R[i][j] = 0.;
			}
		}
		return 1;
	}
	
	cout << "uplo " << uplo << " not understood, upper triangle returned!" << endl;
	return 2;
}

// log determinant of matrix C
double matrix_log_determinant ( const vector< vector<double> > & C ) {
	
	if (C.size() != C[0].size()) {
		throw Exception("matrix_log_determinant: Matrix Not Square!");
	}
	
	long N = C.size();                  // dim of the matrix
	vector< vector<double> > R;
	char UPLO = 'U';                    
	int info = matrix_cholesky(C, R, UPLO);
	if(info <= 0) {
		throw Exception("matrix_log_determinant: Cholesky Factorization Failed!");
	}
	
	double log_det = 0.;
	for (size_t i = 0; i < N; ++i) {
		log_det += log(R[i][i]) * 2.;
	}
	return log_det;
	
}

int matrix_inverse_triangular ( const vector< vector<double> > & R,
                                vector< vector<double> > & X,
                                char uplo) {
	if (R.size() != R[0].size()) {
		throw Exception("matrix_inverse_triangle: Matrix Not Square!");
	}
	if ( !( uplo=='U' || uplo=='u' || uplo=='L' || uplo=='l' ) ) {
		throw Exception("matrix_inverse_triangle: uplo not understood!");
	}
	
	size_t N = R.size();
	X.clear();
	X.resize( N, vector<double>(N, 0.) );
	
	if ( uplo == 'U' || uplo == 'u' ) {
		for (long k = 0; k < N; ++k) {
			for (long i = N-1; i >= 0; --i) {
				X[i][k] = (k==i)?(1.):(0.);
				
				for (long j = i+1; j < N; ++j) {
					X[i][k] -= R[i][j] * X[j][k];
				}
				X[i][k] /= R[i][i];
			}
		}
	}
	
	if ( uplo == 'L' || uplo == 'l' ) {
		for (size_t k = 0; k < N; ++k) {
			for (size_t i = 0; i < N; ++i) {
				X[i][k] = (k==i)?(1.):(0.);
				for (size_t j = 0; j < i; ++j) {
					X[i][k] -= R[i][j] * X[j][k];
				}
				X[i][k] /= R[i][i];
			}
		}
	}
	
	return 1;
}

int matrix_inverse_cholesky ( const vector< vector<double> > & A,
                              vector< vector<double> > & invA) {

	if (A.size() != A[0].size()) {
		cerr << "matrix_inverse_cholesky: Matrix Not Square!" << endl;
		return 0;
	}
	
	vector< vector<double> > R, invR, invRt;
	matrix_cholesky(A, R, 'U');
	matrix_inverse_triangular(R, invR, 'U');
	matrix_transpose(invR, invRt);
	matrix_multiplication(invR, invRt, invA);
	
	return 1;
}

void matrix_display (const vector< vector<double> > & A) {
	
	for (size_t i = 0; i < A.size(); ++i) {
		for (size_t j = 0; j < A[0].size(); ++j) {
			cout << setw(15) << setprecision(10) << A[i][j];
		}
		cout << endl;
	}
}
