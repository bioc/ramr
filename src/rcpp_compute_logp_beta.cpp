#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]

////////////////////////////////////////////////////////////////////////////////

// This is a *modified* version of regularised incomplete beta function
// taken from https://github.com/codeplea/incbeta
// which is distributed under the following license:

/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#define STOP 1.0e-8
#define TINY 1.0e-30

static inline double incbeta (double a,                              /* alpha */
                              double b,                               /* beta */
                              double lbf, /* std::log(std::beta(alpha, beta)) */
                              double x)                                  /* x */
{
  /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
  if (x > (a+1.0)/(a+b+2.0))
    return incbeta(b,a,lbf,1-x);      /*Use the fact that beta is symmetrical.*/
  
  /*Find the log first part before the continued fraction.*/
  const double front = a*std::log(x) + b*std::log(1.0-x) - lbf - std::log(a);
  
  /*Use Lentz's algorithm to evaluate the continued fraction.*/
  double f = 1.0, c = 1.0, d = 0.0;
  
  int i, m;
  for (i = 0; i <= 200; ++i) {
    m = i/2;
    
    double numerator;
    if (i == 0) {
      numerator = 1.0; /*First numerator is 1.0.*/
    } else if (i % 2 == 0) {
      numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
    } else {
      numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
    }
    
    /*Do an iteration of Lentz's algorithm.*/
    d = 1.0 + numerator * d;
    if (std::abs(d) < TINY) d = TINY;
    d = 1.0 / d;
    
    c = 1.0 + numerator / c;
    if (std::abs(c) < TINY) c = TINY;
    
    const double cd = c*d;
    f *= cd;
    
    /*Check for stop.*/
    if (std::abs(1.0-cd) < STOP) {
      return front + std::log(f-1.0);
    }
  }
  
  return NA_REAL; /*Needed more loops, did not converge.*/
}


////////////////////////////////////////////////////////////////////////////////


// To be called after rcpp_fit_beta
//
// Function computes logp values using precomputed alpha, beta and
// log(std::beta), and stores them in 'out' (not transposed anymore)
//
// TODO:
//   [ ] OpenMP
//   [?] maybe skip rows where len[r]==0

// [[Rcpp::export]]
int rcpp_compute_logp_beta (Rcpp::List &data)                                   // List output of rcpp_prepare_data
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_raw> raw((SEXP)data.attr("raw_xptr"));                           // flat vector with raw values
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold output values
  Rcpp::XPtr<T_coef> coef((SEXP)data.attr("coef_xptr"));                        // vector with per-row results of rcpp_fit_beta
  
  // fast direct accessors
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  const auto coef_data = coef->data();
  
  for (size_t c=0; c<ncol; c++) {
    const auto raw_first = raw_data + c*nrow;                                   // first element of c-th column in 'raw'
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    for (size_t r=0; r<nrow; r++) {
      const auto coef_first = coef_data + r*NCOEF;                              // first element of 'coef' array
      out_first[r] = incbeta(coef_first[3], coef_first[4], coef_first[5], raw_first[r]); // Regularized Incomplete Beta Function
    }
  }
  
  return 0;
}





