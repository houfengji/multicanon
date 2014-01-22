#ifndef ENSEMBLEMEAN_H
#define ENSEMBLEMEAN_H

#include <vector>

std::vector<double> ensembleMean (const std::vector< std::vector<double> > &);

std::vector<double> ensembleMean(const std::vector<double> & , const size_t );

double logEnsembleMeanLog (const std::vector<double> &);

std::vector<double> logEnsembleMeanLog(const std::vector<double> & , const size_t );


#endif
