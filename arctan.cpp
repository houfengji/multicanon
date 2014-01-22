#include <cmath>
#include "arctan.h"

// This routine is a customized arc tangent function.
// The return value of this function is in a desired range.
double arctan (const double sine, 
               const double cosine) {
	if (cosine >= 0) {
		if (sine >= 0) {
			return atan(sine / cosine);
		}
		else {
			return atan(sine / cosine) + 2 * M_PI;
		}
	}
	else {
		return atan(sine / cosine) + M_PI;
	} 
}
