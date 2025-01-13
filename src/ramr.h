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
const size_t MINNSMPL = 3;                                                      // minimum number of samples with value within (0;1) (i.e., excluding 0, 1 and NaN)
const size_t NCOEF = 8;                                                         // number of coefficient values to compute per genomic position
// always:         [0]number of zeros, [1]number of ones, [2]median
// for IQR:        [3]Q3, [4]Q1, [5]IQR
// for beta (MoM): [3]mean, [4]variance, [5]alpha, [6]beta, [7]log(std::beta)
// for beta (MLE): [3]sample geometric mean, [4]sample geometric mean based on (1 − X), [5]alpha, [6]beta, [7]log(std::beta)
// for binomial:   [3]log probability of 0, [4]log probability of 1, [5-7]same as beta

// TYPEDEFS //
typedef std::vector<unsigned int> T_chr;                                        // vector of chromosomes
typedef std::vector<unsigned int> T_pos;                                        // vector of genomic positions
typedef std::vector<unsigned int> T_str;                                        // vector of genomic strands
typedef std::vector<double> T_raw;                                              // vector of raw values
typedef std::vector<unsigned int> T_cov;                                        // optional vector of coverages
typedef std::vector<double> T_out;                                              // vector of computed values
typedef std::vector<unsigned int> T_len;                                        // vector of number of columns by row in data
typedef std::vector<double> T_coef;                                             // vector to store coefficients: median, IQR, parameters of fitted distribution

// MACRO //
#define isZero(x) ((x) <= DBL_EPSILON)
#define isOne(x) ((x) >= 1-DBL_EPSILON)
#define notZO(x) (!isZero(x) && !isOne(x))

// OpenMP //
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif



