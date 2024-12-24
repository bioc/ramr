#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]

// To be called after rcpp_compute_*
//
// Function thresholds values in 'out', i.e., assigns NaN if condition is FALSE
//
// TODO:
//   [ ] OpenMP
//   [?] maybe skip rows where len[r]==0

template<bool is_xiqr>
int rcpp_filter_threshold (Rcpp::List &data,                                    // List output of rcpp_prepare_data
                           double thr)                                          // threshold
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: either xIQR or p-values)
  
  // fast direct accessors
  const auto out_data = out->data();
  
  for (size_t c=0; c<ncol; c++) {
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    for (size_t r=0; r<nrow; r++) {
      if (is_xiqr) {                                                            // if xIQR values
        if (std::abs(out_first[r]) < thr) out_first[r] = NA_REAL;               // if less than threshold then make it NaN. Comparisons to NaN is always FALSE
      } else{                                                                   // if p-values
        if (out_first[r] > thr) out_first[r] = NA_REAL;                         // if greater than threshold then make it NaN. Comparisons to NaN is always FALSE
      }
    }
  }
  
  return 0;
}


// [[Rcpp::export]]
int rcpp_filter_threshold_xiqr (Rcpp::List &data, double thr)
{
  return rcpp_filter_threshold<true>(data, thr);
}

// [[Rcpp::export]]
int rcpp_filter_threshold_pval (Rcpp::List &data, double thr)
{
  return rcpp_filter_threshold<false>(data, thr);
}


