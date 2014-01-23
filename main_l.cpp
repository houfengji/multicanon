#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "linearalgebra.h"

int main(void) {
	vector< vector<double> > A(3, vector<double>(3, 0));
	A[0][0] = 9; A[0][1] = 3; A[0][2] = 1;
	A[1][0] = 3; A[1][1] = 7.11212; A[1][2] = 3;
	A[2][0] = 1; A[2][1] = 3; A[2][2] = 5;
	vector< vector<double> > invA;
	matrix_inverse_cholesky(A, invA);
	matrix_display(A);
	matrix_display(invA);
	return 1;
}
