#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]

// To be called after rcpp_get_iqr
//
// Function computes xIQR values using precomputed median and IQR:
//   1) subtracts median from each 'raw' value, divides by IQR
//   2) stores xIQR values in 'out' (not transposed anymore)
//
// TODO:
//   [ ] OpenMP
//   [?] maybe skip rows where len[r]==0

// [[Rcpp::export]]
int rcpp_compute_xiqr (Rcpp::List &data)                                        // List output of rcpp_prepare_data
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_raw> raw((SEXP)data.attr("raw_xptr"));                           // flat vector with raw values
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: xIQR)
  Rcpp::XPtr<T_coef> coef((SEXP)data.attr("coef_xptr"));                        // vector with per-row results of rcpp_get_iqr ([0]median, [1]Q3, [2]Q1, [3]IQR)
  
  // fast direct accessors
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  const auto coef_data = coef->data();
  
  for (size_t c=0; c<ncol; c++) {
    const auto raw_first = raw_data + c*nrow;                                   // first element of c-th column in 'raw'
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    for (size_t r=0; r<nrow; r++) {
      const auto coef_first = coef_data + r*NCOEF;                              // first element of 'coef' array
      out_first[r] = (raw_first[r] - coef_first[0]) / coef_first[3];            // (value-median)/IQR
    }
  }
  
  return 0;
}

