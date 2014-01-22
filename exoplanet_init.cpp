#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "exoplanet_hyperpara.h"
#include "exoplanet_init.h"

using namespace std;

ExoplanetInit::ExoplanetInit(string data_name, size_t num_comp, const ExoplanetHyperpara & hyper) {
	
	bool found = 0;
	
	if (data_name.compare("data1") == 0 && num_comp == 0) {
		found = 1;
		
		init.resize(num_comp * 5 + 2);
		
		init[10] = 2.5;
		init[11] = 6000.;
	}
	
	if (data_name.compare("data1") == 0 && num_comp == 1) {
		found = 1;
		
		init.resize(num_comp * 5 + 2);
		
		init[0] = 100.;
		init[1] = 2*M_PI/100.;
		init[2] = 1.;
		init[3] = 0.05;
		init[4] = 4.;
		init[10] = 2.;
		init[11] = 100.;
	}
	
	if (data_name.compare("data1") == 0 && num_comp == 2) {
		found = 1;
		
		init.resize(num_comp * 5 + 2);
		
		init[0] = 10.;
		init[1] = 2*M_PI/30.;
		init[2] = 2.2;
		init[3] = 0.5;
		init[4] = 5;
		init[5] = 100.;
		init[6] = 2*M_PI/100.;
		init[7] = 1;
		init[8] = 0.05;
		init[9] = 4;
		init[10] = 2.;
		init[11] = 30.;
	}
	
	
	if (found == 0) {
		cerr << "Exoplanet Initialization Failed!" << endl;
		cerr << "Data Name or Model couldn't be found!" << endl;
	}
}

ExoplanetInit::ExoplanetInit(string data_name, size_t num_comp, size_t num_d, const ExoplanetHyperpara & hyper) {
	
	bool found = 0;
	
	if (found == 0) {
		cerr << "Exoplanet Initialization Failed!" << endl;
		cerr << "Data Name or Model couldn't be found!" << endl;
	}
}
