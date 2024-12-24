// FLOW:
//   1) rcpp_prepare_data <- transposes 'raw' values into 'out'
//   2) rcpp_(get|fit)_(iqr|beta|beinf) <- gets median/IQR or fits a
//      distribution; 'coef' stores these values
//   3) rcpp_compute_(xiqr|pbeta|pbeinf) <- computes xIQR values or p-values
//      from 'raw' values and 'coef' values, stores in 'out'
//   4) rcpp_filter_(threshold|combp) <- combines values in regions,
//      either applying a threshold or by comb-p approach
//   5) rcpp_create_granges <- makes GRanges object with AMRs
//

// CONSTS //
const size_t ncoef = 4;                                                         // number of coefficient values to compute per genomic position


// TYPEDEFS //
typedef std::vector<uint32_t> T_chr;                                            // vector of chromosomes
typedef std::vector<uint32_t> T_pos;                                            // vector of genomic positions
typedef std::vector<uint32_t> T_str;                                            // vector of genomic strands
typedef std::vector<double> T_raw;                                              // vector of raw values
typedef std::vector<double> T_out;                                              // vector of computed values
typedef std::vector<uint32_t> T_len;                                            // vector of number of columns by row in data
typedef std::vector<double> T_coef;                                             // vector to store coefficients: median, IQR, parameters of fitted distribution


