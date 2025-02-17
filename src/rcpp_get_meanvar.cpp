#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// Function computes mean and variance, and stores them
// in the vector of coefficients as one of the following:
//   a) {[3] mean(x), [4] var(x)}
//   b) {[3] weighted mean(x), [4] weighted var(x)}
//   c) {[3] geometric mean(x), [4] geometric mean(1-x)}
//   d) {[3] weighted geometric mean(x), [4] weighted geometric mean(1-x)}

// REFS:
//  1) https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
//  2) https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
//  3) https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Variance-defined_weights
//  4) https://en.wikipedia.org/wiki/Inverse-variance_weighting

// TODO:
//   [x] make it ready for 0 and 1 - now it is not aware of them
//   [ ] OpenMP
//   [ ] ...

// MACRO //
#define invDist(x) (1 / (std::abs(q[2] - (x)) + FLT_EPSILON))                   /* weight inversely correlates with distance from the median */
#define sqrtInvDist(x) (std::sqrt(invDist(x)))                                  /* weight inversely correlates with square root of distance from the median */
#define logInvDist(x) (std::log(invDist(x)))                                    /* weight is a negative logarithm of distance from the median */
#define getWeight(x) {                                                         \
  switch (tWeight) {                                                           \
  case 0:                                                  /* equal weights */ \
    break;                                                                     \
  case 1:                                /* weight = 1 / abs(x - median(x)) */ \
    w = invDist(x); break;                                                     \
  case 2:                       /* weight = sqrt ( 1 / abs(x - median(x)) ) */ \
    w = sqrtInvDist(x); break;                                                 \
  case 3:                        /* weight = log ( 1 / abs(x - median(x)) ) */ \
    w = logInvDist(x); break;                                                  \
  }                                                                            \
};


template<int tMean, int tWeight>                                                // templated for 'type of mean' and 'type of weighting'
int rcpp_get_meanvar (Rcpp::List &data)                                         // List output of rcpp_prepare_data
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)

  // containers
  Rcpp::XPtr<T_dbl> out((SEXP)data.attr("out_xptr"));                           // vector with intermediate output values (here: transposed 'raw')
  Rcpp::XPtr<T_int> len((SEXP)data.attr("len_xptr"));                           // lengths of input data rows minus number of NaNs
  Rcpp::XPtr<T_dbl> coef((SEXP)data.attr("coef_xptr"));                         // vector to hold per-row results

  // fast direct accessors
  const auto out_data = out->data();
  const auto len_data = len->data();
  const auto coef_data = coef->data();

  for (size_t r=0; r<nrow; r++) {
    const auto first = out_data + r*ncol;                                       // first element
    const auto q = coef_data + r*NCOEF;                                         // pointer to the first element of 'coef' NCOEF-element array
    const size_t l = len_data[r];                                               // length = ncol - nNaNs
    const size_t lzo = (size_t)(q[0]+q[1]+0.5);                                 // number of 0s and 1s within l
    if (l < ((tMean==0 ? 0 : lzo) + MINNSMPL)) {                                // if not enough values to process (arithmetic mean accepts 0/1)
      std::fill_n(q+3, NCOEF-3, NA_REAL);                                       // estimates are NaN
      continue;                                                                 // skip this row
    }
    double w = 1;                                                               // weight
    double sumweights = 0;                                                      // accumulator of weights
    double sumsquares = 0;                                                      // accumulator of weights^2

    // NB: ARITHMETIC MEAN ALLOWS 0/1, GEOMETRIC MEAN SKIPS 0/1
    if (tMean==0) {                                                             // arithmetic mean and variance
      // mean in q[3]
      q[3] = 0;
      for (size_t i=0; i<l; i++) {
        getWeight(first[i]);
        q[3] += first[i] * w;
        sumweights += w;
      }
      q[3] /= sumweights;

      // unbiased variance in q[4]
      q[4] = 0;
      for (size_t i=0; i<l; i++) {
        getWeight(first[i]);
        q[4] += std::pow(first[i] - q[3], 2) * w;
        sumsquares += std::pow(w, 2);
      }
      q[4] /= sumweights - sumsquares/sumweights;                               // https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights
                                                                                // or as in Hmisc::wtd.var(x, w, normwt=TRUE)
    } else if (tMean==1) {                                                      // geometric means
      // sample geometric mean is in q[3]
      // sample geometric mean based on (1 − X) is in q[4]
      q[3] = 0;
      q[4] = 0;
      for (size_t i=0; i<l; i++) {
        if (notZO(first[i])) {
          getWeight(first[i]);
          q[3] += std::log(first[i]) * w;
          q[4] += std::log(1 - first[i]) * w;
          sumweights += w;
        }
      }
      q[3] = exp(q[3] / sumweights);
      q[4] = exp(q[4] / sumweights);
    }
  }

  return 0;
}

// [[Rcpp::export]]
int rcpp_get_meanvar_ari_equal (Rcpp::List &data)                               // arithmetic mean, equal weights
{
  return rcpp_get_meanvar<0, 0>(data);
}

// [[Rcpp::export]]
int rcpp_get_meanvar_ari_invDist (Rcpp::List &data)                             // arithmetic mean, inverse distance weights
{
  return rcpp_get_meanvar<0, 1>(data);
}

// [[Rcpp::export]]
int rcpp_get_meanvar_ari_sqrtInvDist (Rcpp::List &data)                         // arithmetic mean, sqrt inverse distance weights
{
  return rcpp_get_meanvar<0, 2>(data);
}

// [[Rcpp::export]]
int rcpp_get_meanvar_ari_logInvDist (Rcpp::List &data)                          // arithmetic mean, log inverse distance weights
{
  return rcpp_get_meanvar<0, 3>(data);
}

// [[Rcpp::export]]
int rcpp_get_meanvar_geo_equal (Rcpp::List &data)                               // geometric mean, equal weights
{
  return rcpp_get_meanvar<1, 0>(data);
}

// [[Rcpp::export]]
int rcpp_get_meanvar_geo_invDist (Rcpp::List &data)                             // geometric mean, inverse distance weights
{
  return rcpp_get_meanvar<1, 1>(data);
}

// [[Rcpp::export]]
int rcpp_get_meanvar_geo_sqrtInvDist (Rcpp::List &data)                         // geometric mean, sqrt inverse distance weights
{
  return rcpp_get_meanvar<1, 2>(data);
}

// [[Rcpp::export]]
int rcpp_get_meanvar_geo_logInvDist (Rcpp::List &data)                          // geometric mean, log inverse distance weights
{
  return rcpp_get_meanvar<1, 3>(data);
}







