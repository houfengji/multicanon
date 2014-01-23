#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <vector>

using namespace std;

bool symmetric (const vector< vector<double> > & C);

int matrix_multiplication ( vector< vector<double> > & A, vector< vector<double> > & B, vector< vector<double> > & C);

int matrix_transpose( vector< vector<double> > & A, vector< vector<double> > & At);

int matrix_cholesky( const vector< vector<double> > & A,
                     vector< vector<double> > & R,
                     char uplo);

double matrix_log_determinant (const vector< vector<double> > & C);

int matrix_inverse_triangular ( const vector< vector<double> > & R,
                                vector< vector<double> > & X,
                                char uplo);

int matrix_inverse_cholesky ( const vector< vector<double> > & A,
                              vector< vector<double> > & invA);

void matrix_display (const vector< vector<double> > & A);

#endif
