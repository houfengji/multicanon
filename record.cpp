#include <cmath>
#include <fstream>
#include <iomanip>
#include "model.h"
#include "record.h"

using namespace std;

void record ( const Model & model, std::fstream & out) {
	out << setprecision(16) << model.beta.back() << "   " << model.chain_evi.back() << "   " << exp(model.chain_evi.back()) << "   " << model.chain_C.back() << "   " << model.chain_R.back() << "   " << model.tau.back() << endl;
}
