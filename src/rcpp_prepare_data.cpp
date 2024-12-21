// #include <algorithm>
#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp20)]]

// This function prepares input data for further processing:
//   1) makes a copy of raw methylation values ('raw')
//   2) transposes raw values dropping NaNs ('out') and counting them ('len')
//   3) arranges other vectors used in computations later.
//
// TODO:
//   [ ] more efficient access to S4Vectors with raw values
//   [ ] OpenMP
//   [ ] ...

// <input.ranges> for rcpp_prepare_data must be sorted
// [[Rcpp::export]]
Rcpp::List rcpp_prepare_data (Rcpp::IntegerVector &seqnames,                    // IntegerVector output of as.integer(S4Vectors::runValue(GenomeInfoDb::seqnames(<input.ranges>)))
                              Rcpp::IntegerVector &seqrunlens,                  // IntegerVector output of as.integer(S4Vectors::runLength(GenomeInfoDb::seqnames(<input.ranges>)))
                              Rcpp::IntegerVector &start,                       // IntegerVector output of as.integer(BiocGenerics::start(<input.ranges>))
                              Rcpp::IntegerVector &strand,                      // IntegerVector output of as.integer(BiocGenerics::strand(<input.ranges>))
                              Rcpp::DataFrame &mcols)                           // DataFrame output of as.data.frame(GenomicRanges::mcols(<input.ranges>))
{
#define T_coef std::array<double,4>                                             // array to store median, IQR, parameters of fitted distribution, etc
  
  // consts
  const size_t ncol = mcols.ncol();                                             // number of columns (samples)
  const size_t nrow = mcols.nrow();                                             // number of rows (genomic loci)
  
  // containers
  std::vector<uint32_t>* pos = new std::vector<uint32_t>(start.begin(), start.end());    // genomic positions
  std::vector<uint32_t>* str = new std::vector<uint32_t>(strand.begin(), strand.end());  // genomic strands
  std::vector<double>* raw = new std::vector<double>;                           // flat vector with raw values from &mcols
  std::vector<double>* out = new std::vector<double>;                           // vector to hold intermediate output values (e.g., transposed)
  std::vector<uint32_t>* len = new std::vector<uint32_t>;                       // lengths of &mcols rows minus number of NaNs
  std::vector<T_coef>* coef = new std::vector<T_coef>;                          // vector to hold per-row results (e.g., median, Q1, Q3, parameters of fitted distribution)

  // fill 'raw' with values from &mcols
  raw->reserve(ncol*nrow);                                                      // reserve space as required
  for (size_t c=0; c<ncol; c++)
    raw->insert(raw->end(), ((Rcpp::NumericVector)mcols[c]).begin(), ((Rcpp::NumericVector)mcols[c]).end());
  raw->shrink_to_fit();
  
  // initialize 'len', 'coef', and 'out'
  len->resize(nrow);                                                            // init with 0 - a waste, but no choice
  coef->resize(nrow);                                                           // have a feeling that values are not initialized to 0, but that's ok
  out->resize(ncol*nrow, NA_REAL);                                              // init with NA_REAL - see if it breaks anything further. NB: default might be marginally faster
  
  // fast direct accessors
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  const auto len_data = len->data();
  
  // transpose 'raw' to 'out', skipping NaNs; adjust 'len'
  // should be more computationally efficient and parallelizable
  // have to rewrite this to become cache-friendly, 845 samples seriously suck on Mac
  double *buf  = (double*) malloc(ncol * sizeof(double));                       // buffer to gather values from each column (mcols[r,])
  for (size_t r=0; r<nrow; r++) {
    size_t l = 0;                                                               // number of elements actually copied
    for (size_t c=0; c<ncol; c++)                                               // column by column
      if (!std::isnan(raw_data[r+nrow*c]))                                      // if value is not a NaN
        buf[l++] = raw_data[r+nrow*c];                                          // gather it in the buffer; increase its length
    len_data[r] = l;                                                            // adjust observed length
    std::copy(buf, buf+l, out_data+ncol*r);                                     // copy 'buf' to 'out'
  }
  free(buf);
  
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

devtools::clean_dll()
pkgbuild::compile_dll(debug=FALSE)
devtools::load_all()

devtools::check()
cvg <- covr::package_coverage(type="all")
cvg; covr::zero_coverage(cvg)

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

rcpp_test_nan()

load("~/work/data/ramr/data/GSE51032/GSE51032.data.Rdata")
test.ranges <- geo.ranges

load("~/work/data/ramr/data/REVISION-SIMULATED/REVISION-SIMULATED-5-1000-0.250.data.Rdata")
test.ranges <- simulated.ranges

test.ranges <- multi.ranges
seqnames_ <- as.integer(S4Vectors::runValue(GenomeInfoDb::seqnames(test.ranges)))
seqrunlens_ <- as.integer(S4Vectors::runLength(GenomeInfoDb::seqnames(test.ranges)))
start_ <- as.integer(BiocGenerics::start(test.ranges))
strand_ <- as.integer(BiocGenerics::strand(test.ranges))
mcols_ <- as.data.frame(GenomicRanges::mcols(test.ranges, use.names=FALSE))

microbenchmark::microbenchmark(
  z <- rcpp_prepare_data(seqnames_, seqrunlens_, start_, strand_, mcols_),
  y <- t(mcols_),
  {suppressWarnings(rm(y,z)); gc(); gc()},
times=25)

microbenchmark::microbenchmark(
  y <- as.data.frame(GenomicRanges::mcols(test.ranges)),
  z <- as.numeric(as(GenomicRanges::mcols(test.ranges, use.names=FALSE), "Vector")),
  {suppressWarnings(rm(y,z)); gc(); gc()},
  times=10)

*/
// #############################################################################

