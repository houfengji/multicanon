/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *
 */
 
#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include <cstdlib>
#include <string>
#include <vector>
#include "data.h"
#include "model.h"

class Rosenbrock: public Model {
public:
	Rosenbrock();
	double x1_lower_bound;  // hyper-parameter
	double x1_upper_bound;  // hyper-parameter
	double x2_lower_bound;  // hyper-parameter
	double x2_upper_bound;  // hyper-parameter
	virtual double reparametrize(const std::vector<double> &, std::vector<double> &) const;
	virtual double LnDensity(const std::vector<double> &) ;
	virtual double LnLikelihood(const std::vector<double> &);
	virtual double LnImportance(const std::vector<double> &);
	virtual double LnPosterior(const std::vector<double> &);
	double fxx(const std::vector<double> &);  // 2nd order derivative of posterior
	double fyy(const std::vector<double> &);  // 2nd order derivative of posterior
	double LnPrior(const std::vector<double> &) ;
	double LnG(const std::vector<double> &) ;
	void init(size_t, std::vector< std::vector<double> > &, double) const;
	virtual ~Rosenbrock();
};

double rosenBruteIntegral(Rosenbrock & rosen, const size_t mesh_size_1, const size_t mesh_size_2);

#endif
