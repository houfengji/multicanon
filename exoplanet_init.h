#ifndef EXOPLANET_INIT_H
#define EXOPLANET_INIT_H

#include <string>
#include <vector>
#include "exoplanet_hyperpara.h"

struct ExoplanetInit {
	ExoplanetInit(std::string data_name, size_t num_comp, size_t num_d, const ExoplanetHyperpara & hyper);
	ExoplanetInit(std::string data_name, size_t num_comp, const ExoplanetHyperpara & hyper);
	
	std::vector<double> init;
};


#endif
