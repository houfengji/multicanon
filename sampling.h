#ifndef SAMPLING_H
#define SAMPLING_H

#include "model.h"

                        
size_t sampling (Model & model,      // data and model
                 std::vector< std::vector<double> > & ensemble,
                 const double a);
                 
                  
#endif
