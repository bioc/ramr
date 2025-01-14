#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// Function estimates log probability of zeros and ones and stores them
// in the vector of coefficients as {[3] log(p(0)), [4] log(p(1))}
//
// there are alpha and beta available at this point
// so they can be used to estimate the probability of getting all 0s or all 1s
// but it is easier to use average methylation value 'm'
// and have p(0) = 1 - m and p(1) = m
//
// TODO:
//   [ ] OpenMP
//   [ ] ...

// [[Rcpp::export]]
int rcpp_fit_binom (Rcpp::List &data)                                           // List output of rcpp_prepare_data
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector with intermediate output values (here: transposed 'raw')
  Rcpp::XPtr<T_len> len((SEXP)data.attr("len_xptr"));                           // lengths of input data rows minus number of NaNs
  Rcpp::XPtr<T_coef> coef((SEXP)data.attr("coef_xptr"));                        // vector to hold per-row results
  
  // fast direct accessors
  const auto out_data = out->data();
  const auto len_data = len->data();
  const auto coef_data = coef->data();
  
  for (size_t r=0; r<nrow; r++) {
    const auto q = coef_data + r*NCOEF;                                         // pointer to the first element of 'coef' NCOEF-element array
    const size_t l = len_data[r];                                               // length = ncol - nNaNs
    if (isZero(q[0]+q[1]) || l==0) continue;                                    // if no 0/1 or no data, skip this row
    const auto first = out_data + r*ncol;                                       // first element
    
    // mean in m
    double m = 0;
    for (size_t i=0; i<l; i++)
      m += first[i];
    m /= l;
    
    // log probability of 0 is in q[3]
    q[3] = std::log(1 - m);
          
    // log probability of 1 is in q[4]
    q[4] = std::log(m);
  }
  
  return 0;
}



