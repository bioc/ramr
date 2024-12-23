// Flow:
//   1) rcpp_prepare_data <- transposes 'raw' values into 'out'
//   2) rcpp_(get|fit)_(iqr|beta|beinf) <- gets median/IQR or fits a
//      distribution; 'coef' stores these values
//   3) rcpp_compute_(xiqr|pbeta|pbeinf) <- computes xIQR values or p-values
//      from 'raw' values and 'coef' values, stores in 'out'
//   4) rcpp_filter_(threshold|combp) <- combines values in regions,
//      either applying a threshold or by comb-p approach
//   5) rcpp_create_granges <- makes GRanges object with AMRs
//

// DEFINITIONS //

// vector of genomic position
#define T_pos std::vector<uint32_t>

// vector of genomic strands
#define T_str std::vector<uint32_t>

// vector of raw values
#define T_raw std::vector<double>

// vector of genomic strands
#define T_out std::vector<double>

// vector of number of columns in data
#define T_len std::vector<uint32_t>

// number of coefficient values to compute per genomic position
#define ncoef 4

// vector to store coefficients: median, IQR, parameters of fitted distribution
#define T_coef std::vector<double>


