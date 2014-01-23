#include <cmath>
#include <iostream>
#include "factorial.h"

double factorial(long x) {
	if (x < 0) {
		return -1;
	}
	if (x == 0) {
		return 1;
	}
	double r = 1;
	for (long i = 1; i <= x; ++i) {
		r *= i;
	}
	return r;
}

