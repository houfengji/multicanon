/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *
 */
 
#ifndef EXOPLANETJD_H
#define EXOPLANETJD_H


#include <cstdlib>
#include <string>
#include <vector>
#include "data.h"
#include "exoplanet_hyperpara.h"
#include "model.h"

class ExoplanetJD: public Model {
public:
	ExoplanetJD(Data &, size_t num_companion);
	size_t num_comp;    // number of companions
	size_t num_d;       // number of distinctions
	size_t count_distinction();
	double LnPermutation;
	double LnNorm;      // Overall normalization factor
	ExoplanetHyperpara hyper;
	virtual double reparametrize(const std::vector<double> &, std::vector<double> &) const;
	void reverse_reparametrize(const std::vector<double> &, std::vector<double> &) const;
	virtual double LnDensity(const std::vector<double> &);
	virtual double LnImportance(const std::vector<double> &);
	virtual double LnPosterior(const vector<double> &);
	double LnPrior01(const std::vector<double> &);
	double LnPrior02(const std::vector<double> &);
	virtual double LnLikelihood(const std::vector<double> &);
	double LnLikelihoodConstOne(const std::vector<double> &);
	double LnG(const std::vector<double> &) ;
	bool orbit_cross (const std::vector<double> &);
	bool physical (const std::vector<double> &);
	bool exchange (const std::vector<double> &);
	void init(size_t ens_size, std::vector< std::vector<double> > & ensemble, double ini);
	void init(size_t ens_size, std::vector< std::vector<double> > & ensemble, std::vector<double> & V, double ini);
	virtual ~ExoplanetJD();
private:
	Data data;
};

#endif
