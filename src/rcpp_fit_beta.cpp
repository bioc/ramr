#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// Function estimates parameters of beta distribution and stores them
// in the vector of coefficients as {[5] alpha (p), [6] beta (q), [7] log(std::beta)}
//
// TODO:
//   [ ] make it ready for 0 and 1 - now it is not aware of them
//   [ ] OpenMP
//   [ ] ...

template<int method>
int rcpp_fit_beta (Rcpp::List &data)                                            // List output of rcpp_prepare_data
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
    const auto first = out_data + r*ncol;                                       // first element
    const auto q = coef_data + r*NCOEF;                                         // pointer to the first element of 'coef' NCOEF-element array
    const size_t l = len_data[r];                                               // length = ncol - nNaNs
    const size_t lzo = (size_t)(q[0]+q[1]+0.5);                                 // number of 0s and 1s within l
    if (l < (lzo+MINNSMPL)) {                                                   // if not enough values to process (either excluded by median, or too many NaN/0/1)
      std::fill_n(q+3, NCOEF-3, NA_REAL);                                       // estimates are NaN
      continue;                                                                 // skip this row
    }

    // WHICH ONES OF THESE NEED TO SKIP 0/1?
    if (method==0) {                                                            // method of moments
      // mean in q[3]
      q[3] = 0;
      for (size_t i=0; i<l; i++)
        q[3] += first[i];
      q[3] /= l;

      // variance in q[4]
      q[4] = 0;
      for (size_t i=0; i<l; i++)
        q[4] += std::pow(first[i] - q[3], 2);
      q[4] /= l - 1;

      // alpha (shape parameter p) in q[5]
      q[5] = q[3] * (( (q[3] * (1 - q[3])) / q[4]) - 1);

      // beta (shape parameter q) in q[6]
      q[6] = (1 - q[3]) * (((q[3] * (1 - q[3])) / q[4]) - 1);

    } else if (method==1) {                                                     // approximate MLE
      // https://en.wikipedia.org/wiki/Beta_distribution#Maximum_likelihood

      // sample geometric mean in q[3]
      // sample geometric mean based on (1 − X) in q[4]
      q[3] = 0;
      q[4] = 0;
      for (size_t i=0; i<l; i++) {
        q[3] += std::log(first[i]);
        q[4] += std::log(1 - first[i]);
      }
      q[3] = exp(q[3]/l);
      q[4] = exp(q[4]/l);

      // alpha (shape parameter p) in q[5]
      q[5] = 0.5 + q[3] / ( 2 * (1 - q[3] - q[4]) );

      // beta (shape parameter q) in q[6]
      q[6] = 0.5 + q[4] / ( 2 * (1 - q[3] - q[4]) );

    } else if (method==2) {                                                     // TODO: numerical MLE
      Rcpp::stop("not implemented");
      // check how stats::optim works, maybe look for a C++ solution
      // possibly should use multi-objective optimization of the set of my
      // two equations with digamma:
      // https://scicomp.stackexchange.com/questions/3318/simultaneous-maximization-of-two-functions-without-available-derivatives
    }

    // logarithm of complete beta function in q[4]
    q[7] = std::lgamma(q[5]) + std::lgamma(q[6]) - std::lgamma(q[5] + q[6]);
    // there's absolutely no error handling here...
    // but all numbers are defined and finite, so should be fine?..
  }

  return 0;
}

// [[Rcpp::export]]
int rcpp_fit_beta_mom (Rcpp::List &data)                                        // method of moments
{
  return rcpp_fit_beta<0>(data);
}

// [[Rcpp::export]]
int rcpp_fit_beta_amle (Rcpp::List &data)                                       // approximate MLE
{
  return rcpp_fit_beta<1>(data);
}

// [[Rcpp::export]]
int rcpp_fit_beta_nmle (Rcpp::List &data)                                       // numerical MLE
{
  return rcpp_fit_beta<2>(data);
}


