#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]

// To be called after rcpp_get_iqr
//
// Function computes xIQR values using precomputed median and IQR:
//   1) subtracts median from each 'raw' value, divides by IQR
//   2) stores xIQR values in 'out' (not transposed anymore)
//
// TODO:
//   [x] OpenMP
//   [?] maybe skip rows where len[r]==0

// [[Rcpp::export]]
int rcpp_compute_xiqr (Rcpp::List &data)                                        // List output of rcpp_prepare_data
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)

  // containers
  Rcpp::XPtr<T_dbl> raw((SEXP)data.attr("raw_xptr"));                           // flat vector with raw values
  Rcpp::XPtr<T_dbl> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: xIQR)
  Rcpp::XPtr<T_dbl> coef((SEXP)data.attr("coef_xptr"));                         // vector with per-row results of rcpp_get_iqr ([0]median, [1]Q3, [2]Q1, [3]IQR)
  Rcpp::XPtr<T_int> thr((SEXP)data.attr("thr_xptr"));                           // chunks of rows for multiple threads

  // fast direct accessors
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  const auto coef_data = coef->data();

  // number of chunks/threads
  const size_t nthreads = thr->size() - 1;                                      // 'thr' always starts with 0 and ends with 'nrow'

#pragma omp parallel num_threads(nthreads)
{
  const size_t thr_num = omp_get_thread_num();                                  // thread ID
  const size_t row_from = thr->at(thr_num);                                     // start of row chunk
  const size_t row_to = thr->at(thr_num+1);                                     // end of row chunk

  for (size_t c=0; c<ncol; c++) {
    const auto raw_first = raw_data + c*nrow;                                   // first element of c-th column in 'raw'
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    for (size_t r=row_from; r<row_to; r++) {
      const auto coef_first = coef_data + r*NCOEF;                              // first element of 'coef' array
      out_first[r] = (raw_first[r] - coef_first[2]) / coef_first[5];            // (value-median)/IQR
    }
  }
}

  return 0;
}

