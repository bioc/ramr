#include <algorithm>
#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp20)]]

// This function computes xIQR values:
//   1) computes Q1, median and Q3
//   2) subtracts median from each 'raw' value, abs it, divides by IQR
//
// TODO:
//   [ ] OpenMP
//   [ ] ...

// [[Rcpp::export]]
int rcpp_compute_iqr (Rcpp::List &data)                                         // List output of rcpp_prepare_data
{
#define T_coef std::array<double,4>                                             // array to store median, IQR, parameters of fitted distribution, etc
  
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<std::vector<double>> raw((SEXP)data.attr("raw_xptr"));             // flat vector with raw values
  Rcpp::XPtr<std::vector<double>> out((SEXP)data.attr("out_xptr"));             // vector to hold intermediate output values (here: transposed)
  Rcpp::XPtr<std::vector<uint32_t>> len((SEXP)data.attr("len_xptr"));           // lengths of input data rows minus number of NaNs
  Rcpp::XPtr<std::vector<T_coef>> coef((SEXP)data.attr("coef_xptr"));           // vector to hold per-row results (here: [0]median, [1]IQR)
  
  // fast direct accessors
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  const auto len_data = len->data();
  const auto coef_data = coef->data();
  
  for (size_t r=0; r<nrow; r++) {
    const auto first = out_data + r*ncol;                                       // first element
    const size_t l = len_data[r];                                               // length = ncol - NaNs
    std::nth_element(first, first + l*3/4, first+l);
    if ((l/2)&1) {}
  }
  
  return (int)nrow * (int)ncol;
}

