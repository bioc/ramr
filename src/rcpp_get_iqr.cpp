#include <algorithm>
#include <vector>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]

// Function computes xIQR values using R's default quantile function (type 7):
//   1) computes Q1, median and Q3
//   2) stores its output in the vector of coefficients as {IQR, Q1, median, Q3}
//   2) subtracts median from each 'raw' value, abs it, divides by IQR
//
// TODO:
//   [ ] OpenMP
//   [ ] ...

// [[Rcpp::export]]
int rcpp_get_iqr (Rcpp::List &data)                                             // List output of rcpp_prepare_data
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_raw> raw((SEXP)data.attr("raw_xptr"));                           // flat vector with raw values
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: transposed)
  Rcpp::XPtr<T_len> len((SEXP)data.attr("len_xptr"));                           // lengths of input data rows minus number of NaNs
  Rcpp::XPtr<T_coef> coef((SEXP)data.attr("coef_xptr"));                        // vector to hold per-row results (here: [0]IQR, [1]Q1, [2]median, [3]Q3)
  
  // fast direct accessors
  const auto out_data = out->data();
  const auto len_data = len->data();
  const auto coef_data = coef->data();
  
  // arrays for quantile calculations
  const double p[2] = {0.75, 0.25};                                             // probabilities
  double g[2];                                                                  // gamma = g
  size_t j[2];                                                                  // index j
  
  for (size_t r=0; r<nrow; r++) {
    const auto first = out_data + r*ncol;                                       // first element
    const size_t l = len_data[r];                                               // length = ncol - nNaNs
    if (l==0) continue;                                                         // skip row if no values to process
    const auto q = coef_data + r*NCOEF + 1;                                     // pointer to the second element of 'coef' NCOEF-element array
    
    // for Type 7 quantile function (R's default):
    //   m = 1 - p
    //   j = floor(np + m)
    //   gamma = g = np + m - j
    //   Q(p) = (1 - g)Xj + gX(j+1)
    // where p=0.25 for Q1 and p=0.75 for Q3,
    // and Xj and X(j+1) are j-th and j+1-th statistics (elements)
    
    // for Q1 and Q3 quantiles, corner case is when ((l-1)&3)==0 (i.e., 21, 25, 27, ..., 4n+1),
    // then gamma = 0 and only one nth_element() call per Q is therefore required
    // can test it by checking if g<0.1
    
    // for median, corner case is when (l&1)==0 (i.e., 21, 23, 25, ..., 2n+1),
    // then gamma = 0 and only one nth_element() call per Q is therefore required
    
    //    | q1		      |  q2		     |  q3		
    //----|-----|-------|-----|------|-----|------
    // n	|  j	| g	    |  j	| g	   | j	 | g	
    //----|-----|-------|-----|------|-----|------
    // 20	|  5	| 0,75	| 10  | 0,5	 | 15	 | 0,25	
    // 21	|  6	| 0	    | 11	| 0	   | 16	 | 0	
    // 22	|  6	| 0,25	| 11	| 0,5	 | 16	 | 0,75	
    // 23	|  6	| 0,5	  | 12	| 0	   | 17	 | 0,5	
    // 24	|  6	| 0,75	| 12	| 0,5	 | 18	 | 0,25	
    // 25	|  7	| 0	    | 13	| 0	   | 19	 | 0	
    // 26	|  7	| 0,25	| 13	| 0,5	 | 19	 | 0,75	
    // 27	|  7	| 0,5	  | 14	| 0	   | 20	 | 0,5	
    // 28	|  7	| 0,75	| 14	| 0,5	 | 21	 | 0,25	
    // 29	|  8	| 0	    | 15	| 0	   | 22	 | 0	
    // 30	|  8	| 0,25	| 15	| 0,5	 | 22	 | 0,75	
    
    // compute indexes and coefficients for quantiles
    for (size_t i=0; i<2; i++) {
      g[i] = (double)l * p[i] + 1 - p[i];                                       // gamma = g = np + m = np + 1 - p
      j[i] = (size_t)g[i];                                                      // j = floor(g)
      g[i] -= j[i];                                                             // g = g - j
      j[i]--;                                                                   // because must be 0-based, not 1-
    }
    
    // Q3:
    if (g[0]>0.1) {                                                             // ((l-1)&3)!=0, i.e., not the case of g=0
      std::nth_element(first, first+j[0]+1, first+l);                           // get (j+1)-th element
      q[0] = first[j[0]+1] * g[0];                                              // (j+1)-th times gamma
      std::nth_element(first, first+j[0], first+j[0]+1);                        // get j-th element
      q[0] += first[j[0]] * (1-g[0]);                                           // plus j-th times 1-gamma
    } else {
      std::nth_element(first, first+j[0], first+l);                             // get j-th element
      q[0] = first[j[0]];                                                       // which is a Q3
    }
    
    // Q1:
    if (g[1]>0.1) {                                                             // ((l-1)&3)!=0, i.e., not the case of g=0
      std::nth_element(first, first+j[1]+1, first+j[0]);
      q[1] = first[j[1]+1] * g[1];
      std::nth_element(first, first+j[1], first+j[1]+1);
      q[1] += first[j[1]] * (1-g[1]);
    } else {
      std::nth_element(first, first+j[1], first+j[0]);
      q[1] = first[j[1]];
    }
    
    // IQR:
    q[2] = q[0] - q[1];
  }
  
  return 0;
}

