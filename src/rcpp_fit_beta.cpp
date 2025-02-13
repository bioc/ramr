#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// Function estimates parameters of (optionally, weighted) beta distribution
// and stores them in the vector of coefficients as
// {[5] alpha (p), [6] beta (q), [7] log(std::beta)}
//
// TODO:
//   [ ] numerical MLE
//   [ ] OpenMP
//   [ ] ...

template<int method>
int rcpp_fit_beta (Rcpp::List &data)                                            // List output of rcpp_prepare_data
{
  // consts
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)

  // containers
  Rcpp::XPtr<T_coef> coef((SEXP)data.attr("coef_xptr"));                        // vector to hold per-row results

  // fast direct accessors
  const auto coef_data = coef->data();

  for (size_t r=0; r<nrow; r++) {
    const auto q = coef_data + r*NCOEF;                                         // pointer to the first element of 'coef' NCOEF-element array
    if (std::isnan(q[3]))                                                       // if means were not computed (not enough values to process)
      continue;                                                                 // skip this row

    // NB: MOM ALLOWS 0/1, *MLE SKIP 0/1
    if (method==0) {                                                            // method of moments based on the unbiased estimator of variance
      // after rcpp_get_meanvar():
      //   mean in q[3]
      //   variance in q[4]

      // alpha (shape parameter p) is in q[5]
      q[5] = q[3] * (( (q[3] * (1 - q[3])) / q[4]) - 1);

      // beta (shape parameter q) is in q[6]
      q[6] = (1 - q[3]) * (((q[3] * (1 - q[3])) / q[4]) - 1);

    } else if (method==1) {                                                     // approximate MLE
      // https://en.wikipedia.org/wiki/Beta_distribution#Maximum_likelihood
      // after rcpp_get_meanvar():
      //   sample geometric mean is in q[3]
      //   sample geometric mean based on (1 − X) is in q[4]

      // alpha (shape parameter p) is in q[5]
      q[5] = 0.5 + q[3] / ( 2 * (1 - q[3] - q[4]) );

      // beta (shape parameter q) is in q[6]
      q[6] = 0.5 + q[4] / ( 2 * (1 - q[3] - q[4]) );

    } else if (method==2) {                                                     // TODO: numerical MLE
      Rcpp::stop("not implemented");
      // check how stats::optim works: grep -RIn "optim" src/appl/* src/include/*
      // maybe look for a C++ solution
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


