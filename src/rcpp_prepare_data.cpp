#include <algorithm>
#include <ranges>
#include <vector>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(BH)]]

// 


// [[Rcpp::export]]
double rcpp_test ()
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
  std::ranges::nth_element(v, v.begin() + v.size()/2);
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "na-aware nth element:\n";
  std::ranges::nth_element(v, v.begin() + v.size()/2, customLess);
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "default sort:\n";
  std::ranges::sort(v);
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "na-aware sort:\n";
  std::ranges::sort(v, customLess);
  for (int i=0; i<n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  return(0);
}

// <input.ranges> for rcpp_prepare_data must be sorted
// [[Rcpp::export]]
Rcpp::List rcpp_prepare_data (Rcpp::IntegerVector &seqnames,                    // IntegerVector output of as.integer(S4Vectors::runValue(GenomeInfoDb::seqnames(<input.ranges>)))
                              Rcpp::IntegerVector &seqrunlens,                  // IntegerVector output of as.integer(S4Vectors::runLength(GenomeInfoDb::seqnames(<input.ranges>)))
                              Rcpp::IntegerVector &start,                       // IntegerVector output of as.integer(BiocGenerics::start(<input.ranges>))
                              Rcpp::IntegerVector &strand,                      // IntegerVector output of as.integer(BiocGenerics::strand(<input.ranges>))
                              Rcpp::DataFrame &mcols)                           // DataFrame output of GenomicRanges::mcols(<input.ranges>)
{
#define T_coef std::array<double,4>                                             // array to store median, IQR, parameters of fitted distribution, etc
  
  // consts
  const size_t ncol = mcols.ncol();                                             // number of columns (samples)
  const size_t nrow = mcols.nrow();                                             // number of rows (genomic loci)
  
  // containers
  std::vector<uint32_t>* pos = new std::vector<uint32_t>(start.begin(), start.end());    // genomic positions
  std::vector<uint32_t>* str = new std::vector<uint32_t>(strand.begin(), strand.end());  // genomic strands
  std::vector<double>* raw = new std::vector<double>;                           // flat vector with raw values from &mcols
  std::vector<double>* out = new std::vector<double>;                           // vector to hold intermediate output values (e.g., transposed, sorted)
  std::vector<uint32_t>* len = new std::vector<uint32_t>;                       // lengths of &mcols rows minus number of NaNs
  std::vector<T_coef>* coef = new std::vector<T_coef>;                          // vector to hold per-row results (e.g., median, Q1, Q3, parameters of fitted distribution)

  // fill 'raw' with values from &mcols
  raw->reserve(ncol*nrow);                                                      // reserve space as required
  for (size_t c=0; c<ncol; c++)
    raw->insert(raw->end(), ((Rcpp::NumericVector)mcols[c]).begin(), ((Rcpp::NumericVector)mcols[c]).end());
  raw->shrink_to_fit();
  
  // initialize 'len', 'coef', and 'out'
  len->resize(nrow, ncol);                                                      // all values = number of columns (samples)
  coef->resize(nrow);                                                           // have a feeling that values are not initialized to 0
  out->resize(ncol*nrow);                                                       // init with 0 - a waste, but no choice
  
  // transpose 'raw' to 'out', counting NaNs, subtracting them from 'len' vector
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  const auto len_data = len->data();
  for (size_t c=0; c<ncol; c++) {
    for (size_t r=0; r<nrow; r++) {
      const double v = raw_data[nrow*c+r];
      out_data[ncol*r+c] = v;
      if (std::isnan(v)) len_data[r]--;
    }
  }
  
  // sort 'out' putting NaNs on the right
  struct {
    bool operator()(double a, double b) const {return std::isnan(b) || (a<b);}  // NaN last
  } nanLess;
  for (size_t r=0; r<nrow; r++) {
    std::sort(out_data+r*ncol, out_data+(r+1)*ncol, nanLess);
  }
  
  // 
  
  
  // wrap and return the results
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("ncol") = ncol,                                                 // number of columns (samples)
    Rcpp::Named("nrow") = nrow,                                                 // number of rows (genomic loci)
    Rcpp::Named("seqnames") = seqnames,                                         // integer IDs of seqnames
    Rcpp::Named("seqrunlens") = seqrunlens                                      // running lengths of seqnames
  );
  
  // pointers to containers
  Rcpp::XPtr<std::vector<uint32_t>> pos_xptr(pos, true);
  Rcpp::XPtr<std::vector<uint32_t>> str_xptr(str, true);
  Rcpp::XPtr<std::vector<double>> raw_xptr(raw, true);
  Rcpp::XPtr<std::vector<double>> out_xptr(out, true);
  Rcpp::XPtr<std::vector<uint32_t>> len_xptr(len, true);
  Rcpp::XPtr<std::vector<T_coef>> coef_xptr(coef, true);
  res.attr("pos_xptr") = pos_xptr;
  res.attr("str_xptr") = str_xptr;
  res.attr("raw_xptr") = raw_xptr;
  res.attr("out_xptr") = out_xptr;
  res.attr("len_xptr") = len_xptr;
  res.attr("coef_xptr") = coef_xptr;
  
  return(res);
}



// #############################################################################
// test code and sourcing don't work on OS X
/*** R
options(width=140)
setwd("~/work/packages/ramr/")
devtools::document()
devtools::load_all()

devtools::check()

# library(data.table)
library(GenomicRanges)
data(ramr)
sn <- seqnames(ramr.data)
multi.ranges <- unlist(as(lapply(levels(sn), function (chr) {
  S4Vectors::runValue(sn) <- factor(chr, levels=levels(sn))
  gr <- ramr.data
  seqnames(gr) <- sn
  return(gr)
}), "GRangesList"))

S4Vectors::runLength(seqnames(multi.ranges))
S4Vectors::runValue(seqnames(multi.ranges))

rcpp_test()

seqnames_ <- as.integer(S4Vectors::runValue(GenomeInfoDb::seqnames(multi.ranges)))
seqrunlens_ <- as.integer(S4Vectors::runLength(GenomeInfoDb::seqnames(multi.ranges)))
start_ <- as.integer(BiocGenerics::start(multi.ranges))
strand_ <- as.integer(BiocGenerics::strand(multi.ranges))
mcols_ <- GenomicRanges::mcols(multi.ranges)

microbenchmark::microbenchmark(
  z <- rcpp_prepare_data(seqnames_, seqrunlens_, start_, strand_, mcols_),
times=100)


*/
// #############################################################################

