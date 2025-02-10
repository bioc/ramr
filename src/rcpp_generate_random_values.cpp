#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// Function generates random values using estimated probability (inferred from
// observed counts) of zeros (in [0]) and ones (in [1]), and estimated
// parameters of beta distribution (in [5] and [6]).
// Outputs 1D vector of doubles (by column).
//
// TODO:
//   [ ] OpenMP, but find a way to seed properly
//   [ ] ...

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_generate_random_values (Rcpp::List &data,              // List output of rcpp_prepare_data
                                                 const size_t ncol)             // number of samples
{
  // consts
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)

  // containers
  Rcpp::XPtr<T_len> len((SEXP)data.attr("len_xptr"));                           // lengths of input data rows minus number of NaNs
  Rcpp::XPtr<T_coef> coef((SEXP)data.attr("coef_xptr"));                        // vector to hold per-row results

  // result
  std::vector<double> res(nrow*ncol, NA_REAL);                                  // instantiate with NA

  // fast direct accessors
  const auto len_data = len->data();
  const auto coef_data = coef->data();
  const auto res_data = res.data();

  for (size_t r=0; r<nrow; r++) {
    const auto q = coef_data + r*NCOEF;                                         // pointer to the first element of 'coef' NCOEF-element array
    if (len_data[r]==0) continue;                                               // skip this row if empty
    const double p0 = q[0] / len_data[r];                                       // probability of 0
    const double p1 = q[1] / (len_data[r] - q[0]);                              // probability of 1 among non-0
    const bool has0 = !isZero(q[0]);                                            // are there 0s?
    const bool has1 = !isZero(q[1]);                                            // are there 1s?
    const bool hasB = !std::isnan(q[5]);                                        // are beta estimates valid?

    for (size_t c=0; c<ncol; c++) {                                             // sample by sample
      if (has0 && R::rbinom(1, p0)>0.5) {                                       // if probability of 0s is >0 and it's a 0
        res_data[r+c*nrow] = 0;
      } else if (has1 && R::rbinom(1, p1)>0.5) {                                // else if probability of 1s is >0 and it's a 1
        res_data[r+c*nrow] = 1;
      } else if (hasB) {                                                        // if neither zero nor one and estimates are valid
        res_data[r+c*nrow] = R::rbeta(q[5], q[6]);                              // random beta
      }
    }
  }

  Rcpp::NumericVector res_matrix = Rcpp::wrap(res);                             // wrap it
  Rcpp::IntegerVector dim = {(int)nrow, (int)ncol};
  res_matrix.attr("dim") = dim;                                                 // set dimensions
  return res_matrix;
}



