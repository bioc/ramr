#include <vector>
#include <functional>
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


// [[Rcpp::export]]
int rcpp_filter_threshold_xiqr (Rcpp::List &data,                               // List output of rcpp_prepare_data
                                double thr)                                     // threshold
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: xIQR)
  
  // fast direct accessors
  const auto out_data = out->data();
  
  for (size_t c=0; c<ncol; c++) {
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    for (size_t r=0; r<nrow; r++) {
      if (std::abs(out_first[r]) < thr) out_first[r] = NA_REAL;                 // if TRUE then make it NaN. Comparisons to NaN is always FALSE
    }
  }
  
  return 0;
}

// [[Rcpp::export]]
int rcpp_filter_threshold_pval (Rcpp::List &data,                               // List output of rcpp_prepare_data
                                double thr)                                     // threshold
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: xIQR)
  
  // fast direct accessors
  const auto out_data = out->data();
  
  for (size_t c=0; c<ncol; c++) {
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    for (size_t r=0; r<nrow; r++) {
      if (out_first[r] > thr) out_first[r] = NA_REAL;                           // if TRUE then make it NaN. Comparisons to NaN is always FALSE
    }
  }
  
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// templated version is slow...
//
// 
// template<typename T_Tfm, typename T_Cmp>
// int rcpp_filter_threshold (Rcpp::List &data,                                    // List output of rcpp_prepare_data
//                            double thr,                                          // threshold
//                            T_Tfm Tfm,                                           // transform function: identity (for p-values) or abs (for xIQR)
//                            T_Cmp Cmp)                                           // comparison function: less (for xIQR) or greater (for p-values)
// {
//   // consts
//   const size_t ncol = data["ncol"];                                             // number of columns (samples)
//   const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
//   
//   // containers
//   Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: xIQR)
//   
//   // fast direct accessors
//   const auto out_data = out->data();
//   
//   for (size_t c=0; c<ncol; c++) {
//     const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
//     for (size_t r=0; r<nrow; r++) {
//       if (Cmp(Tfm(out_first[r]), thr)) out_first[r] = NA_REAL;                  // if TRUE then make it NaN. Comparisons to NaN is always FALSE
//     }
//   }
//   
//   return 0;
// }
// 
// // [[Rcpp::export]]
// int rcpp_filter_threshold_xiqr (Rcpp::List &data, double thr)
// {
//   return rcpp_filter_threshold(data, thr, std::abs<double>, std::less<double>());
// }
// 
// // [[Rcpp::export]]
// int rcpp_filter_threshold_pval (Rcpp::List &data, double thr)
// {
//   return rcpp_filter_threshold(data, thr, std::abs<double>, std::greater<double>());  // have to make std::identity<double>() work
// }



