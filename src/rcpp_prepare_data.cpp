// #include <algorithm>
#include <vector>
#include <array>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]

// This function prepares input data for further processing:
//   1) makes a copy of raw methylation values ('raw')
//   2) transposes raw values dropping NaNs (to 'out')
//      and counting 0s, 1s and other valid values ('coef', 'len')
//   3) arranges other vectors used in computations later.
//
// TODO:
//   [ ] more efficient access to S4Vectors with raw values
//   [ ] OpenMP
//   [ ] ...

// <input.ranges> for rcpp_prepare_data must be sorted
// [[Rcpp::export]]
Rcpp::List rcpp_prepare_data (Rcpp::IntegerVector &seqnames,                    // IntegerVector (factor) output of S4Vectors::runValue(GenomeInfoDb::seqnames(<input.ranges>))
                              Rcpp::IntegerVector &seqrunlens,                  // IntegerVector output of S4Vectors::runLength(GenomeInfoDb::seqnames(<input.ranges>))
                              Rcpp::IntegerVector &start,                       // IntegerVector output of BiocGenerics::start(<input.ranges>)
                              Rcpp::IntegerVector &strand,                      // IntegerVector (factor) output of as.factor(BiocGenerics::strand(<input.ranges>))
                              Rcpp::DataFrame &mcols,                           // DataFrame output of as.data.frame(GenomicRanges::mcols(<input.ranges>))
                              double exclude_lower,                             // lower bound of range to exclude
                              double exclude_upper)                             // upper bound of range to exclude
{
  // consts
  const size_t ncol = mcols.ncol();                                             // number of columns (samples)
  const size_t nrow = mcols.nrow();                                             // number of rows (genomic loci)

  // containers
  T_chr* chr = new T_chr;                                                       // chromosomes
  T_pos* pos = new T_pos(start.begin(), start.end());                           // genomic positions
  T_str* str = new T_str(strand.begin(), strand.end());                         // genomic strands
  T_raw* raw = new T_raw;                                                       // flat vector with raw values from &mcols
  T_out* out = new T_out;                                                       // vector to hold intermediate output values (e.g., transposed)
  T_len* len = new T_len;                                                       // lengths of &mcols rows minus number of NaNs
  T_coef* coef = new T_coef;                                                    // vector to hold per-row results (e.g., median, Q1, Q3, parameters of fitted distribution)

  // fill 'chr' vector with seqname ids
  chr->reserve(nrow);                                                           // reserve space as required
  for (size_t i=0; i<seqnames.size(); i++)
    chr->resize(chr->size()+seqrunlens[i], seqnames[i]);
  chr->shrink_to_fit();

  // fill 'raw' with values from &mcols
  raw->reserve(ncol*nrow);                                                      // reserve space as required
  for (size_t c=0; c<ncol; c++)
    raw->insert(raw->end(), ((Rcpp::NumericVector)mcols[c]).begin(), ((Rcpp::NumericVector)mcols[c]).end());
  raw->shrink_to_fit();

  // initialize 'len', 'coef', and 'out'
  len->resize(nrow);                                                            // init with 0
  coef->resize(nrow*NCOEF);                                                     // nrow times NCOEF to store them all continuously (all 0)
  out->resize(ncol*nrow, NA_REAL);                                              // init with NA_REAL - see if it breaks anything further. NB: default might be marginally faster

  // fast direct accessors
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  const auto len_data = len->data();
  const auto coef_data = coef->data();

  // transpose 'raw' to 'out', counting 0/1, skipping NaNs; adjust 'len'
  // should be more computationally efficient and parallelizable
  // have to rewrite this to become cache-friendly, 845 samples seriously suck on Mac
  double *buf  = (double*) malloc(ncol * sizeof(double));                       // buffer to gather values from each column (mcols[r,])
  for (size_t r=0; r<nrow; r++) {
    const auto q = coef_data + r*NCOEF;                                         // pointer to coef NCOEF-element array
    size_t l = 0;                                                               // number of elements actually copied
    for (size_t c=0; c<ncol; c++) {                                             // column by column
      const auto raw_value = raw_data[r+nrow*c];                                // value to compare/transpose
      if (!std::isnan(raw_value)) {                                             // if is not NaN
        q[0] += isZero(raw_value);                                              // is it a 0?
        q[1] += isOne(raw_value);                                               // is it a 1?
        buf[l++] = raw_value;                                                   // gather it in the buffer; increase its length
      }
    }

    // median
    if (l>0) {                                                                  // if there are values in the buffer
      const size_t hl = l/2;                                                    // half length
      std::nth_element(buf, buf+hl, buf+l);                                     // order up to l/2-th
      q[3] = buf[hl];                                                           // median for odd l
      if ((l&1)==0) {                                                           // if l is even
        std::nth_element(buf, buf+hl-1, buf+hl);                                // order up to l/2-1-th
        q[3] = (q[3] + buf[hl-1])/2;                                            // median for even l
      }
      if (q[3]<exclude_lower || q[3]>exclude_upper){                            // if median is less that exclude_lower or greater than exclude_upper
        std::copy(buf, buf+l, out_data+ncol*r);                                 // copy 'buf' to 'out'
        len_data[r] = l;                                                        // adjust observed length, because otherwise it's 0 and we won't use this row in further analyses
      }
    }
  }
  free(buf);

  // wrap and return the results
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("ncol") = ncol,                                                 // number of columns (samples)
    Rcpp::Named("nrow") = nrow,                                                 // number of rows (genomic loci)
    Rcpp::Named("seqnames") = seqnames,                                         // integer IDs of seqnames
    Rcpp::Named("seqrunlens") = seqrunlens,                                     // running lengths of seqnames
    Rcpp::Named("samples") = mcols.names()                                      // sample names
  );
  res.attr("strandlevels") = strand.attr("levels");                             // strand levels

  // pointers to containers
  Rcpp::XPtr<T_chr> chr_xptr(chr, true);
  Rcpp::XPtr<T_pos> pos_xptr(pos, true);
  Rcpp::XPtr<T_str> str_xptr(str, true);
  Rcpp::XPtr<T_raw> raw_xptr(raw, true);
  Rcpp::XPtr<T_out> out_xptr(out, true);
  Rcpp::XPtr<T_len> len_xptr(len, true);
  Rcpp::XPtr<T_coef> coef_xptr(coef, true);
  res.attr("chr_xptr") = chr_xptr;
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

*/
// #############################################################################

