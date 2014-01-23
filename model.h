/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  Feb 22, 2013
 *
 */
 
#ifndef MODEL_H
#define MODEL_H

#include <cstdlib>
#include <string>
#include <vector>

class Model {
public:
	Model(size_t dimension);         // the dimension of the model's paramter space
	Model(void);
	size_t dim;                      // dimension of the model's parameter space
	std::vector<double> best_fit;    // This vector is used to keep track of the best fit parameters
	std::string model_name;          // name of the model, for output file names
	std::string time_label;          // the time when an object is made, used to label output files
	std::vector<double> beta;        // inverse of temperature
	std::vector<double> chain_evi;   // chain of evidence (log of)
	std::vector<double> chain_C;     // chain of the factor controling error bar
	std::vector<double> chain_R;     // chain of the actual relative error bar
	std::vector<double> tau;         // chain of auto-correlation time
	double bic_evi;                  // evidence value from BIC
	virtual double reparametrize(const std::vector<double> &, std::vector<double> &) const = 0;
	virtual double LnDensity(const std::vector<double> &) = 0;
	virtual double LnImportance(const std::vector<double> &) = 0;
	virtual double LnLikelihood(const std::vector<double> &) = 0;
	virtual double LnPosterior(const std::vector<double> &) = 0;
	std::vector< std::vector<double> > C;  // covariance
	double log_det_C;
	std::vector< std::vector<double> > H;  // inverse of covariance
	std::vector<double> m;       // mean
	char Uplo;
	void displayMeanCov(void);
	void displayHessian(void);
	virtual ~Model() {}
private:
};

#endif
