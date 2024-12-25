#include <vector>
#include <list>
#include <Rcpp.h>
#include "ramr.h"

// [[Rcpp::plugins(cpp20)]]

// To be called after rcpp_filter_*
//
// Function assembles a list structure with data needed to create a
// GenomicRanges object with AMRs.
//
// TODO:
//   [ ] OpenMP
//   [ ] ...
//

// [[Rcpp::export]]
Rcpp::List rcpp_create_granges (Rcpp::List &data,                               // List output of rcpp_prepare_data
                                size_t window,                                  // 
                                bool ignore_strand,                             // 
                                size_t min_ncpg,                                // 
                                size_t min_width)                               // 
{
  // consts
  const size_t ncol = data["ncol"];                                             // number of columns (samples)
  const size_t nrow = data["nrow"];                                             // number of rows (genomic loci)
  
  // containers
  Rcpp::XPtr<T_chr> chr((SEXP)data.attr("chr_xptr"));                           // chromosomes (1-based)
  Rcpp::XPtr<T_pos> pos((SEXP)data.attr("pos_xptr"));                           // genomic positions (1-based)
  Rcpp::XPtr<T_str> str((SEXP)data.attr("str_xptr"));                           // genomic strands (1=="+", 2=="-", 3=="*")
  Rcpp::XPtr<T_raw> raw((SEXP)data.attr("raw_xptr"));                           // flat vector with raw values
  Rcpp::XPtr<T_out> out((SEXP)data.attr("out_xptr"));                           // vector to hold intermediate output values (here: xIQR or p-values or NaN for the ones to skip)
  
  // fast direct accessors
  const auto chr_data = chr->data();
  const auto pos_data = pos->data();
  const auto str_data = str->data();
  const auto raw_data = raw->data();
  const auto out_data = out->data();
  
  // output containers for AMRs
  // have to be careful with them when writing from multiple threads
  T_chr res_chr;                                                                // chromosomes
  T_pos res_start;                                                              // genomic start
  T_pos res_end;                                                                // genomic end
  T_str res_strand;                                                             // genomic strand
  std::list<T_pos> res_revmap;                                                  // revmap
  T_pos res_ncpg;                                                               // number of CpGs
  T_pos res_sample;                                                             // integer sample id
  T_raw res_dbeta;                                                              // average 'raw' minus 'median' (beta)
  T_out res_aggr;                                                               // average 'out' (mean for xIQR, geometric mean for p-values), or comb-p combined p
  
  // macros
#define spit_amr {             /* save AMR when enough CpGs and wide enough */ \
  if ((amr[s].revmap.size()>=min_ncpg) &&          /* if ncpg>=min_ncpg and */ \
      (amr[s].end-amr[s].start+1>=min_width)) {         /* width>=min_width */ \
    res_chr.push_back(amr[s].chr);                            /* chromosome */ \
    res_start.push_back(amr[s].start);                             /* start */ \
    res_end.push_back(amr[s].end);                                   /* end */ \
    res_strand.push_back(s+1);                                    /* strand */ \
    res_revmap.push_back(amr[s].revmap);                          /* revmap */ \
    res_ncpg.push_back(amr[s].revmap.size());                       /* ncpg */ \
    amr[s].dbeta /= amr[s].revmap.size();                  /* average dbeta */ \
    res_dbeta.push_back(amr[s].dbeta);                     /* average dbeta */ \
    amr[s].aggr /= amr[s].revmap.size();         /* aggregated 'out' values */ \
    res_aggr.push_back(amr[s].aggr);             /* aggregated 'out' values */ \
  }                                                                            \
  amr[s].revmap.clear();                                    /* clear revmap */ \
  amr[s].open = false;                                         /* close AMR */ \
};

  // cycle through genomic position
  for (size_t c=0; c<ncol; c++) {
    const auto raw_first = raw_data + c*nrow;                                   // first element of c-th column in 'raw'
    const auto out_first = out_data + c*nrow;                                   // first element of c-th column in 'out'
    
    // three structures to hold AMR data for every strand
    struct {                                                                    // for every strand:
      bool open = false;                                                        //   AMR range was opened
      size_t chr;                                                               //   AMR range chromosome
      size_t start;                                                             //   AMR range start
      size_t end;                                                               //   AMR range end
      T_pos revmap;                                                             //   vector to hold revmap
      double dbeta;                                                             //   dbeta
      double aggr;                                                              //   aggregated 'out' values
    } amr[3];
    
    size_t s;                                                                   // strand holder
    
    for (size_t r=0; r<nrow; r++) {
      if (std::isnan(out_first[r])) continue;                                   // next if NaN
      s = str_data[r] - 1;                                                      // strand of current position, 0-based
      if (amr[s].open) {                                                        // if there's an open AMR on this strand
        const size_t d = pos_data[r] - amr[s].end;                              //   distance from previous base
        if ((d<=window) && (amr[s].chr==chr_data[r])) {                         //   if within the window and the same chromosome
          amr[s].end = pos_data[r];                                             //     new end
          amr[s].revmap.push_back(r);                                           //     add element to revmap
          amr[s].dbeta += raw_first[r];                                         //     add 'raw'
          amr[s].aggr += out_first[r];                                          //     add 'out'
        } else {                                                                //   if outside the window or another chromosome
          spit_amr;                                                             //     save existing
        }
      }
      if (!amr[s].open) {                                                       // if we are at the beginning of new AMR (also because we just saved some)
        amr[s].open = true;                                                     //   open it
        amr[s].chr = chr_data[r];                                               //   record chromosome
        amr[s].start = pos_data[r];                                             //   start
        amr[s].end = pos_data[r];                                               //   end = start
        amr[s].revmap.push_back(r);                                             //   add element to revmap
        amr[s].dbeta = raw_first[r];                                            //   first 'raw'
        amr[s].aggr = out_first[r];                                             //   first 'out'
      }
    }
    
    // save last AMR(s) after cycling through all genomic position
    for (s=0; s<3; s++)
      if (amr[s].open) spit_amr;
    
    // add (the same) sample id for all sample AMRs
    res_sample.resize(res_chr.size(), c);
  }
  
  
  
  // wrap and return the results
  Rcpp::IntegerVector col_chr = Rcpp::wrap(res_chr);                            // making seqnames a factor
  col_chr.attr("class") = "factor";
  col_chr.attr("levels") = ((Rcpp::IntegerVector)(data["seqnames"])).attr("levels");
  
  Rcpp::IntegerVector col_strand = Rcpp::wrap(res_strand);                      // making strand a factor
  col_strand.attr("class") = "factor";
  col_strand.attr("levels") = data.attr("strandlevels");
  
  Rcpp::IntegerVector col_start = Rcpp::wrap(res_start);                        // int start
  Rcpp::IntegerVector col_end = Rcpp::wrap(res_end);                            // int start
  Rcpp::IntegerVector col_start = Rcpp::wrap(res_start);                        // int start
  Rcpp::IntegerVector col_ncpg = Rcpp::wrap(res_ncpg);                          // int start
  
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("seqnames") = col_chr,                                          // chromosomes
    Rcpp::Named("start") = col_start,                                           // genomic start
    Rcpp::Named("end") = col_end,                                               // genomic end
    Rcpp::Named("strand") = col_strand,                                         // genomic strand
    Rcpp::Named("revmap") = res_revmap,                                         // revmap
    Rcpp::Named("ncpg") = col_ncpg,                                             // number of CpGs
    Rcpp::Named("sample") = res_sample,                                         // integer sample id
    Rcpp::Named("dbeta") = res_dbeta,                                           // average 'raw' minus 'median' (beta)
    Rcpp::Named("xIQR") = res_aggr                                              // average 'out' (mean for xIQR, geometric mean for p-values), or comb-p combined p
  );
  
  return(res);
}



