#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]

// To be called after rcpp_compute_*
//
// Function thresholds values in 'out', i.e., assigns NaN if condition is FALSE
//
// TODO:
//   [x] OpenMP
//   [?] maybe skip rows where len[r]==0

template<bool is_xiqr>
int rcpp_filter_threshold (Rcpp::List &data,                                    // List output of rcpp_prepare_data
                           double thrshld)                                      // threshold
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)

  // containers
  Rcpp::XPtr<T_dbl> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: either xIQR or p-values)
  Rcpp::XPtr<T_int> thr((SEXP)data.attr("thr_xptr"));                           // chunks of rows for multiple threads

  // fast direct accessors
  const auto out_data = out->data();

  // number of chunks/threads
  const size_t nthreads = thr->size() - 1;                                      // 'thr' always starts with 0 and ends with 'nrow'

#pragma omp parallel num_threads(nthreads)
{
  const size_t thr_num = omp_get_thread_num();                                  // thread ID
  const size_t row_from = thr->at(thr_num);                                     // start of row chunk
  const size_t row_to = thr->at(thr_num+1);                                     // end of row chunk

  for (size_t c=0; c<ncol; c++) {
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    for (size_t r=row_from; r<row_to; r++) {
      if (is_xiqr) {                                                            // if xIQR values
        if (std::abs(out_first[r]) < thrshld) out_first[r] = NA_REAL;           // if less than threshold then make it NaN. Comparisons to NaN is always FALSE
      } else{                                                                   // if p-values
        if (out_first[r] > thrshld) out_first[r] = NA_REAL;                     // if greater than threshold then make it NaN. Comparisons to NaN is always FALSE
      }
    }
  }
}

  return 0;
}


// [[Rcpp::export]]
int rcpp_filter_threshold_xiqr (Rcpp::List &data, double thrshld)
{
  return rcpp_filter_threshold<true>(data, thrshld);
}

// [[Rcpp::export]]
int rcpp_filter_threshold_logp (Rcpp::List &data, double thrshld)
{
  return rcpp_filter_threshold<false>(data, thrshld);
}


