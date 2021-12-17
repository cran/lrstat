#include <Rcpp.h>
using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

void set_seed(int seed);

NumericVector stl_sort(NumericVector x);

IntegerVector findInterval2(NumericVector x, NumericVector breaks);
 
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol);

#endif // __UTILITIES__
