#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
int rcpp_test_omp()
{
  int res = -1; // not available
#if defined(_OPENMP)
  res = omp_get_max_threads();
#endif
  return res;
}
