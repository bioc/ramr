#include <algorithm>
#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// 
  
  
// [[Rcpp::export]]
double rcpp_test_nan ()
{
  std::vector<double> v = {5, 6, 4, 3, 2, 6, 7, 9, 3, 1, 2, 4, NA_REAL, R_NaN, R_PosInf, R_NegInf, 4};
  
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





