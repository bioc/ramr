#include <algorithm>
#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// 

// [[Rcpp::export]]
std::vector<double> rcpp_extract_out (Rcpp::List &data)
{
  Rcpp::XPtr<std::vector<double>> out((SEXP)data.attr("out_xptr"));
  std::vector<double> res(out->begin(), out->end());
  return res;
}
  
// [[Rcpp::export]]
double rcpp_test_nan ()
{
  // set.seed(1); paste(sample(23), collapse=", ")
  std::vector<double> v = {4, 7, 1, 2, 11, 14, 21, 5, 16, 10, 6, 18, 22, 9, 15, 12, 17, 19, NA_REAL, R_NaN, R_PosInf, R_NegInf, 8, 13, 20, 3, 23};
  
  int n = v.size();
  for (int i = 0; i < n; ++i) {
    if (std::isnan(v[i]))
      Rprintf("v[%i] is std::isnan.\n", i);
    if (Rcpp::NumericVector::is_na(v[i]))
      Rprintf("v[%i] is NA.\n", i);
    if (Rcpp::traits::is_nan<REALSXP>(v[i]))
      Rprintf("v[%i] is NaN.\n", i);
    if (Rcpp::traits::is_infinite<REALSXP>(v[i]))
      Rprintf("v[%i] is Inf or -Inf.\n", i);
  }
  
  struct {
    bool operator()(double a, double b) const {return std::isnan(b) || (a<b);}  // NaN last
  } customLess;
  
  Rcpp::Rcout << "unsorted:\n";
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "default nth element:\n";
  std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "na-aware nth element:\n";
  std::nth_element(v.begin(), v.begin() + v.size()/2, v.end(), customLess);
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "default sort:\n";
  std::sort(v.begin(), v.end());
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "na-aware sort:\n";
  std::sort(v.begin(), v.end(), customLess);
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  return(0);
}

#include <boost/math/statistics/univariate_statistics.hpp>
// [[Rcpp::export]]
double rcpp_test_med_boost (std::vector<double> v)
{
  const auto first = v.data();
  const size_t l = v.size();
  double m = boost::math::statistics::median(first, first+l);
  return(m);
}
// [[Rcpp::export]]
double rcpp_test_iqr_boost (std::vector<double> v)
{
  const auto first = v.data();
  const size_t l = v.size();
  double i = boost::math::statistics::interquartile_range(first, first+l);
  return(i);
}
// [[Rcpp::export]]
double rcpp_test_med (std::vector<double> v)
{
  const auto first = v.data();
  const size_t l = v.size();
  const size_t n = l/2;
  std::nth_element(first, first+n, first+l);
  double m = first[n];
  if ((l&1)==0) {
    std::nth_element(first, first+n-1, first+l);
    m += first[n-1];
    m /= 2;
  }
  return(m);
}
// [[Rcpp::export]]
std::vector<double> rcpp_test_iqr_type7 (std::vector<double> v)
{
  const auto first = v.data();
  const size_t l = v.size();
  
  // for Type 7 quantile function (R's default):
  //   m = 1 - p
  //   j = floor(np + m)
  //   gamma = g = np + m - j
  //   Q(p) = (1 - g)Xj + gX(j+1)
  // where p=0.25 for Q1 and p=0.75 for Q3,
  // and Xj and X(j+1) are j-th and j+1-th statistics (elements)
  
  // corner case is when (l-1)&3==0 (i.e., 21, 25, 27, ..., 4n+1),
  // then gamma = 0 and only one nth_element() call per Q is therefore required
  // can test it by checking if g<0.1
  
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
  
  double p[4] = {0.00, 0.25, 0.50, 0.75};
  double g[4];
  size_t j[4];
  double q[4];
  
  for (size_t i=1; i<4; i++) {
    g[i] = (double)l * p[i] + 1 - p[i];
    j[i] = (size_t)g[i];
    g[i] -= j[i];
    j[i]--;                                                                     // because must be 0-based, not 1-
  }
  
  // Q3:
  if (g[3]>0.1) {                                                               // ((l-1)&3)!=0, i.e., not the case of g=0
    std::nth_element(first, first+j[3]+1, first+l);
    q[3] = first[j[3]+1] * g[3];
    std::nth_element(first, first+j[3], first+j[3]+1);
    q[3] += first[j[3]] * (1-g[3]);
  } else {
    std::nth_element(first, first+j[3], first+l);
    q[3] = first[j[3]];
  }
  
  // Q2:
  if (g[2]>0.1) {                                                               // (l&1)!=0, i.e., not the case of g=0
    std::nth_element(first, first+j[2]+1, first+j[3]);
    q[2] = first[j[2]+1] * g[2];
    std::nth_element(first, first+j[2], first+j[2]+1);
    q[2] += first[j[2]] * (1-g[2]);
  } else {
    std::nth_element(first, first+j[2], first+j[3]);
    q[2] = first[j[2]];
  }

  // Q1:
  if (g[1]>0.1) {                                                               // ((l-1)&3)!=0, i.e., not the case of g=0
    std::nth_element(first, first+j[1]+1, first+j[2]);
    q[1] = first[j[1]+1] * g[1];
    std::nth_element(first, first+j[1], first+j[1]+1);
    q[1] += first[j[1]] * (1-g[1]);
  } else {
    std::nth_element(first, first+j[1], first+j[2]);
    q[1] = first[j[1]];
  }
  
  std::vector<double> res(q+1, q+4);
  
  return(res);
}
/*
set.seed(1)
vs <- lapply(20:100, sample)
vs <- lapply(1:100, sample, size=10000, replace=TRUE)
vs <- lapply(20:100, rnorm)
all.equal(
  sapply(vs, median), sapply(vs, rcpp_test_med_boost),
  tolerance=1e-15
)
all.equal(
  sapply(vs, IQR), sapply(vs, rcpp_test_iqr_boost),
  tolerance=1e-15
)
all.equal(
  sapply(vs, IQR), sapply(vs, quantile, prob=0.75, names=F)-sapply(vs, quantile, prob=0.25, names=F),
  tolerance=1e-15
)
all.equal(
  sapply(vs, median), sapply(vs, rcpp_test_med),
  tolerance=1e-15
)
all.equal(
  sapply(vs, function (x) {unname(c(quantile(x, probs=0.25), median(x), quantile(x, probs=0.75)))}),
  sapply(vs, rcpp_test_iqr_type7),
  tolerance=1e-30
)
*/  




  


// [[Rcpp::depends(S4Vectors)]]
#include <S4Vectors_interface.h>
// [[Rcpp::export]]
double rcpp_test_s4v (SEXP x)
{
  // Rcpp::Rcout << get_List_elementType(x) << std::endl;
  return(0);
}



////////////////////////////////////////////////////////////////////////////////
// SNIPPETS ////////////////////////////////////////////////////////////////////

// ******** rcpp_prepare_data *********************************************** //
// // transpose 'raw' to 'out', counting NaNs, subtracting them from 'len' vector
// for (size_t c=0; c<ncol; c++) {
//   for (size_t r=0; r<nrow; r++) {
//     const double v = raw_data[nrow*c+r];
//     out_data[ncol*r+c] = v;
//     if (std::isnan(v)) len_data[r]--;                                         // maybe it's better to count NaNs after sorting. Or don't copy them...
//   }
// }
// 
// // sort 'out' putting NaNs on the right
// struct {
//   bool operator()(double a, double b) const {return std::isnan(b) || (a<b);}  // NaN last
// } nanLess;
// for (size_t r=0; r<nrow; r++) {
//   std::sort(out_data+r*ncol, out_data+(r+1)*ncol, nanLess);
// }

// std::sort(buf, buf+l);                                                      // sort 'buf' - eventually might go for several calls of nth_element()
// std::nth_element(buf, buf+l/2, buf+l);





