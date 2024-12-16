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
  // for (int i = 0; i < n; ++i) {
  //   if (Rcpp::NumericVector::is_na(v[i]))
  //     Rprintf("v[%i] is NA.\n", i);
  //   if (Rcpp::traits::is_nan<REALSXP>(v[i]))
  //     Rprintf("v[%i] is NaN.\n", i);
  //   if (Rcpp::traits::is_infinite<REALSXP>(v[i]))
  //     Rprintf("v[%i] is Inf or -Inf.\n", i);
  // }
  
  for (int i = 0; i < n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  std::ranges::nth_element(v, v.begin() + v.size()/2);
  for (int i = 0; i < n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  std::ranges::sort(v);
  for (int i = 0; i < n; ++i)
    Rcpp::Rcout << v[i] << " ";
  Rcpp::Rcout << "\n";
  
  return(0);
}

// [[Rcpp::export]]
Rcpp::List rcpp_prepare_data (std::string fn,                                    // input: a name of (optionally bgzipped and/or indexed) FASTA file
                              int nthreads)                                      // HTSlib threads, >0 for multiple
{
  // containers
  std::vector<uint64_t> rid;                                                    // numeric ids of reference sequences
  std::vector<std::string> rname;                                               // names of reference sequences
  std::vector<uint64_t> rlen;                                                   // lengths of reference sequences
  std::vector<std::string>* rseq = new std::vector<std::string>;                // reference sequences themselves
  
  // // file IO
  // faidx_t *faidx = fai_load(fn.c_str());                                        // FASTA index
  // if (!faidx) Rcpp::stop("Unable to open FASTA index for reading");             // fall back if error
  // 
  // hts_tpool *tpool = NULL;                                                      // thread pool works correctly
  // if (nthreads>0) {
  //   tpool = hts_tpool_init(nthreads);
  //   fai_thread_pool(faidx, tpool, 0);                                           // new API for FAI thread pool
  //   // bgzf_thread_pool(*(BGZF **)faidx, tpool, 0);                             // old, dirty hack: conversion of faidx_t to BGZF because it's first in the struct
  // }
  // 
  // // vars
  // int flen = 0;                                                                 // fetched sequence length
  // 
  // // fetch sequences
  // for (size_t i=0; i<(unsigned int)faidx_nseq(faidx); i++) {
  //   rid.push_back(i);
  //   const char *name = faidx_iseq(faidx, i);
  //   int64_t length = faidx_seq_len(faidx, name);
  //   
  //   rname.push_back(name);
  //   rlen.push_back(length);
  //   
  //   char *sequence = faidx_fetch_seq(faidx, name, 0, length-1, &flen);          // try fetch sequence
  //   if (length!=flen) Rcpp::stop("Corrupted FASTA index. Delete and try again");// if fetched bytes differ from expected
  //   for (size_t j=0; j<(unsigned int)length; j++)
  //     sequence[j]=acgnt_filter_table[(unsigned char)sequence[j]];               // replace extended IUPAC with N
  //   
  //   rseq->emplace_back((const char*) sequence, length);                         // emplace sequence
  //   free(sequence);                                                             // free allocated
  // }
  // 
  // fai_destroy(faidx);                                                           // free allocated
  // if (tpool) hts_tpool_destroy(tpool);                                          // free thread pool
  
  // wrap and return the results
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("rid") = rid,                                                   // numeric ids of reference sequences
    Rcpp::Named("rname") = rname,                                               // names of reference sequences
    Rcpp::Named("rlen") = rlen                                                  // lengths of reference sequences
  );
  
  Rcpp::XPtr<std::vector<std::string>> rseq_xptr(rseq, true);
  res.attr("rseq_xptr") = rseq_xptr;                                            // external pointer to sequences
  
  return(res);
}


