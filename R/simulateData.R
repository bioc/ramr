#' Template-based simulation of methylation data sets
#'
#' @description
#' `simulateData` generates aberration-free methylation data using an
#' experimental data set as a template, and further introduces methylation
#' aberrations if `GRanges` object containing a set of aberrantly methylated
#' regions was provided. The output can be used to evaluate performance
#' of algorithms for search of differentially (DMR) or aberrantly (AMR)
#' methylated regions.
#'
#' @details
#' For every genomic location in the template data set (`GRanges` object with
#' genomic locations and corresponding methylation beta values included as
#' metadata) `simulateData` does the following:
#' \itemize{
#'   \item estimates parameters of beta distribution
#'   \item in the same input, calculates frequencies of zero and one values
#'   (endpoints; whenever present)
#'   \item uses estimated parameters of beta distribution and probabilities
#'   (observed frequencies) of \{0;1\} values to generate `nsamples` random
#'   values by means of `stats::rbeta` function (for beta values) and/or
#'   `stats::rbinom` function (for \{0;1\} endpoint values, according to their
#'   frequencies and therefore probabilities).
#' }
#' This results in "smoothed" data set that has biologically relevant
#' distribution of methylation values at every genomic location,
#' but does not contain methylation
#' aberrations. If the `amr.ranges` parameter points to a `GRanges` object with
#' aberrations, every AMR is then introduced into the "smoothed" data set as
#' following: if mean methylation beta value for AMR region across all samples
#' in the "smoothed" data set is above (below) 0.5 then all beta values for the
#' sample defined by the `sample` metadata column are decreased (increased) by
#' the absolute value specified in the `dbeta` metadata column. Resulting
#' data sets with (or without) AMR together with the `amr.ranges` set of true
#' positive aberrations can be used as test data set to evaluate performance
#' of algorithms for search of differentially (DMR) or aberrantly (AMR)
#' methylated regions.
#'
#' @note
#' NA values within metadata columns of `template.ranges` are silently dropped
#' in all computations. NA values will also not appear in the result of
#' this function, unless parameters of beta distribution and/or probabilities
#' of zeros or ones cannot be estimated (e.g., due to too many NA values in
#' `template.ranges` metadata).
#'
#' @param template.ranges A `GRanges` object with genomic locations and
#' corresponding methylation beta values included as metadata (the same object
#' must be supplied to this and to the `simulateAMR` functions).
#' @param nsamples A single integer >= 1 indicating the number of samples to
#' generate.
#' @param amr.ranges A `GRanges` object with genomic locations of methylation
#' aberrations (epimutations). If `NULL` (the default), no aberrations is
#' introduced, and function will return "smoothed" data set. If supplied,
#' `GRanges` object must contain the following metadata columns:
#' \itemize{
#'   \item `revmap` -- integer list of `template.ranges` genomic locations that
#'   are included in this AMR region
#'   \item `sample` -- an identifier of a sample to which corresponding AMR
#'   belongs. Must be among the supplied or auto generated `sample.names`
#'   \item `dbeta` -- absolute deviation to be introduced. Must be numeric
#'   within the closed interval [0,1] or NA. When NA - the resulting beta value
#'   for the corresponding genomic position will also be NA
#' }
#' Such an object can be obtained using \code{\link{simulateAMR}} method or
#' manually.
#' @param sample.names A character vector with sample names. If `NULL` (the
#' default), sample names will be auto generated. When specified, the length
#' of the `sample.names` vector must be equal to the value of `nsamples`.
#' @param compute A single string for the distribution to fit to the data.
#' Currently accepts "beta+binom" (the default) only, which stands for
#' endpoint-inflated beta distribution. See Details section
#' and \code{\link{getAMR}} method description for additional explanations.
#' @param compute.estimate A single string for the method of parameter
#' estimation of beta distribution. The default ("mom") stands for the method
#' of moments based on the unbiased estimator of variance and includes \{0;1\}
#' endpoints in calculation of moments (mean, unbiased variance).
#' Other options are "amle" (approximation of maximum likelihood estimation)
#' and "nmle" (numeric maximum likelihood estimation) - both ignore \{0;1\}
#' endpoints in calculations. More details on these methods are given in
#' \code{\link{getAMR}} method description.
#' @param compute.weights A single string for the method to compute optional
#' sample weights that are used during estimation of beta distribution
#' parameters. If default ("equal"), all weights are equal. Otherwise, weight of
#' a value equals to a natural logarithm of inverse absolute distance of this
#' value to the sample median ("logInvDist"), a square root of inverse absolute
#' distance of this value to the sample median ("sqrtInvDist"), or an inverse
#' absolute distance of this value to the sample median ("invDist"). Using
#' weighted parameter estimation allows to increase sensitivity of outlier
#' detection. More details on weighted parameter estimation are given in
#' \code{\link{getAMR}} method description.
#' @param ncores A single integer >= 1 for the number of processes for parallel
#' computation. By default (NULL), function will use half of available cores.
#' Results of this function are identical (reproducible)
#' \strong{only when random seed is set} by `set.seed` function
#' \strong{and a single core is used}.
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return The output is a `GRanges` object with genomic ranges that are equal
#' to the genomic ranges of the provided template and metadata columns
#' containing generated methylation beta values for `nsamples` samples. If
#' `amr.ranges` object was supplied, then randomly generated beta values will be
#' modified accordingly.
#' @seealso \code{\link{simulateAMR}} for the generation of random methylation
#' aberrations, \code{\link{getAMR}} for identification of AMRs,
#' \code{\link{plotAMR}} for plotting AMRs, \code{\link{getUniverse}}
#' for info on enrichment analysis, and `ramr` vignettes for the description of
#' usage and sample data.
#' @examples
#'   data(ramr)
#'   amrs <-
#'     simulateAMR(ramr.data, nsamples=10, regions.per.sample=3,
#'                 samples.per.region=1, min.cpgs=5, merge.window=1000)
#'   noise <-
#'     simulateAMR(ramr.data, nsamples=10, regions.per.sample=20,
#'                 exclude.ranges=amrs, min.cpgs=1, max.cpgs=1, merge.window=1)
#'   noisy.data <-
#'     simulateData(template.ranges=ramr.data, nsamples=10, amr.ranges=c(amrs,noise))
#'   plotAMR(data.ranges=noisy.data, amr.ranges=amrs[1])
#' @export
simulateData <- function (template.ranges,
                          nsamples,
                          amr.ranges=NULL,
                          sample.names=NULL,
                          compute="beta+binom",
                          compute.estimate=c("mom", "amle", "nmle"),
                          compute.weights=c("equal", "logInvDist", "sqrtInvDist", "invDist"),
                          ncores=NULL,
                          verbose=TRUE)
{
  if (!methods::is(template.ranges,"GRanges"))
    stop("'template.ranges' must be a GRanges object")
  if (!is.null(sample.names) & length(sample.names)!=nsamples)
    stop("'sample.names' length must be equal to 'nsamples'")

  if (is.null(sample.names))
    sample.names <- sprintf(paste0("sample%0", nchar(as.character(nsamples)), "i"), seq_len(nsamples))
  if (!is.null(amr.ranges)) {
    if(!methods::is(amr.ranges,"GRanges"))
      stop("'amr.ranges' must be a GRanges object")
    if ( is.null(amr.ranges$revmap) | !all(unlist(amr.ranges$revmap) %in% seq_along(template.ranges)) )
      stop("Malformed 'amr.ranges' object: 'revmap' field is missing or contains out-of-range values")
    if ( is.null(amr.ranges$sample) | !all(amr.ranges$sample %in% sample.names) )
      stop("Malformed 'amr.ranges' object: 'sample' field is missing or its elements are absent from 'sample.names'")
    if ( !is.numeric(amr.ranges$dbeta) | !all(stats::na.omit(amr.ranges$dbeta)>=0) | !all(stats::na.omit(amr.ranges$dbeta)<=1) )
      stop("Malformed 'amr.ranges' object: 'dbeta' field is missing or is outside the closed interval [0,1]")
  }

  template.mcols <- GenomicRanges::mcols(template.ranges)
  template.samples <- colnames(template.mcols)
  compute.estimate <- match.arg(compute.estimate)
  compute.weights <- match.arg(compute.weights)

  #####################################################################################

  .data <- .preprocessData(
    data.ranges=template.ranges, data.samples=template.samples,
    data.coverage=data.frame(), transform="identity",
    exclude.range=c(2,0), ncores=ncores, verbose=verbose
  )

  random.betas <- .getRandomValues(
    data.list=.data, estimate=compute.estimate, weights=compute.weights,
    nsamples=nsamples, sample.names=sample.names, verbose=verbose
  )

  if (!is.null(amr.ranges)) {
    random.betas <- .addEpimutations(
      random.data=random.betas, amr.ranges=amr.ranges, verbose=verbose
    )
  }

  data.ranges <- GenomicRanges::granges(template.ranges)
  GenomicRanges::mcols(data.ranges) <- random.betas

  return(data.ranges)
}
