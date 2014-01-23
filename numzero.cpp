#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace std;

size_t numZero(const vector<double> & X) {
	size_t N = 0;
	for (vector<double>::const_iterator iter = X.begin(); iter != X.end(); ++iter) {
		if (*iter == 0) ++N;
	}
	return N;
}

size_t numZeroLog(const vector<double> & logX) {
	double maxLogX = *max_element(logX.begin(), logX.end());
	size_t N = 0;
	for (vector<double>::const_iterator iter = logX.begin(); iter != logX.end(); ++iter) {
		if ( *iter == -numeric_limits<double>::infinity() ) ++N;
	}
	return N;
}
