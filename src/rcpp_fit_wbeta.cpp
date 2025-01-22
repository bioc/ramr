// Function estimates parameters of *weighted* beta distribution,
// where weights inversely correlate with variance, and stores them
// in the vector of coefficients as {[5] alpha (p), [6] beta (q), [7] log(std::beta)}

// REFS:
//  1) https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
//  2) https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
//  3) https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Variance-defined_weights
//  4) https://en.wikipedia.org/wiki/Inverse-variance_weighting
