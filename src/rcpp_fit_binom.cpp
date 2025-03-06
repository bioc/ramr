#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]

// Function estimates log probability of zeros and ones and stores them
// in the vector of coefficients as {[3] log(p(0)), [4] log(p(1))}
//
// there are alpha and beta available at this point
// so they can be used to estimate the probability of getting all 0s or all 1s
// but it is easier to use average methylation value 'm'
// and have p(0) = 1 - m and p(1) = m
//
// TODO:
//   [x] OpenMP
//   [ ] ...

// [[Rcpp::export]]
int rcpp_fit_binom (Rcpp::List &data)                                           // List output of rcpp_prepare_data
{
  // containers
  Rcpp::XPtr<T_int> len((SEXP)data.attr("len_xptr"));                           // lengths of input data rows minus number of NaNs
  Rcpp::XPtr<T_dbl> coef((SEXP)data.attr("coef_xptr"));                         // vector to hold per-row results
  Rcpp::XPtr<T_int> thr((SEXP)data.attr("thr_xptr"));                           // chunks of rows for multiple threads

  // fast direct accessors
  const auto len_data = len->data();
  const auto coef_data = coef->data();

  // number of chunks/threads
  const size_t nthreads = thr->size() - 1;                                      // 'thr' always starts with 0 and ends with 'nrow'

#pragma omp parallel num_threads(nthreads)
{
  const size_t thr_num = omp_get_thread_num();                                  // thread ID
  const size_t row_from = thr->at(thr_num);                                     // start of row chunk
  const size_t row_to = thr->at(thr_num+1);                                     // end of row chunk

  for (size_t r=row_from; r<row_to; r++) {
    const auto q = coef_data + r*NCOEF;                                         // pointer to the first element of 'coef' NCOEF-element array
    if (isZero(q[0]+q[1]) || len_data[r]==0) continue;                          // if no 0/1 or no data, skip this row

    // mean in q[3] after 'rcpp_get_meanvar'
    const double m = q[3];

    // log probability of 0 goes to q[3]
    q[3] = std::log(1 - m);

    // log probability of 1 goes to q[4]
    q[4] = std::log(m);
  }
}

  return 0;
}



