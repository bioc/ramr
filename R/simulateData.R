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
#' genomic locations and corresponding beta values included as metadata)
#' `simulateData` estimates the parameters of beta distribution by means of
#' `EnvStats::ebeta` function, and then uses estimated parameters to generate
#' `nsamples` random beta values by means of `stats::rbeta` function. This
#' results in "smoothed" data set that has biologically relevant distribution
#' of beta values at every genomic location, but does not contain methylation
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
#' @param template.ranges A `GRanges` object with genomic locations and
#' corresponding beta values included as metadata (same object
#' must be supplied to this and to the `simulateAMR` functions).
#' @param nsamples A single integer >= 1 indicating the number of samples to
#' generate.
#' @param amr.ranges A `GRanges` object with genomic locations of (rare)
#' methylation aberrations. If `NULL` (the default), no aberrations is
#' introduced, and function will return "smoothed" data set. If supplied,
#' `GRanges` object must contain the following metadata columns:
#' \itemize{
#'   \item `revmap` -- integer list of `template.ranges` genomic locations that
#'   are included in this AMR region
#'   \item `sample` -- an identifier of a sample to which corresponding AMR
#'   belongs. Must be among the supplied or auto generated `sample.names`
#'   \item `dbeta` -- absolute deviation to be introduced. Must be numeric
#'   within the range c(0,1) or NA. When NA - the resulting beta value for
#'   the corresponding genomic position will also be NA
#' }
#' Such an object can be obtained using \code{\link{simulateAMR}} method or
#' manually.
#' @param sample.names A character vector with sample names. If `NULL` (the
#' default), sample names will be computed as
#' `paste0("sample", seq_len(nsamples))`. When specified, the length of the
#' `sample.names` vector must be equal to the value of `nsamples`.
#' @param min.beta A single numeric within the range c(0,1). All beta values
#' in the generated data set below this value will be assigned this value.
#' The default: 0.001.
#' @param max.beta A single numeric within the range c(0,1). All beta values
#' in the generated data set above this value will be assigned this value.
#' The default: 0.999.
#' @param cores A single integer >= 1. Number of processes for parallel
#' computation (the default: all but one cores). Results of parallel processing
#' are fully reproducible when the same seed is used (thanks to doRNG).
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
#'     simulateData(ramr.data, nsamples=10, amr.ranges=c(amrs,noise), cores=2)
#'   plotAMR(noisy.data, amr.ranges=amrs[1])
#' @export
simulateData <- function (template.ranges,
                          nsamples,
                          amr.ranges=NULL,
                          sample.names=NULL,
                          compute="beta",
                          compute.estimate=c("mom", "amle", "nmle"),
                          compute.weights=c("equal", "invDist", "sqrtInvDist", "logInvDist"),
                          ncores=NULL,
                          verbose=TRUE)
{
  if (!methods::is(template.ranges,"GRanges"))
    stop("'template.ranges' must be a GRanges object")
  if (!is.null(sample.names) & length(sample.names)!=nsamples)
    stop("'sample.names' length must be equal to 'nsamples'")

  if (is.null(sample.names))
    sample.names <- paste0("sample", seq_len(nsamples))
  if (!is.null(amr.ranges)) {
    if(!methods::is(amr.ranges,"GRanges"))
      stop("'amr.ranges' must be a GRanges object")
    if ( is.null(amr.ranges$revmap) | !all(unlist(amr.ranges$revmap) %in% seq_along(template.ranges)) )
      stop("Malformed 'amr.ranges' object: 'revmap' field is missing or contains out-of-range values")
    if ( is.null(amr.ranges$sample) | !all(amr.ranges$sample %in% sample.names) )
      stop("Malformed 'amr.ranges' object: 'sample' field is missing or its elements are absent from 'sample.names'")
    if ( !is.numeric(amr.ranges$dbeta) | !all(stats::na.omit(amr.ranges$dbeta)>=0) | !all(stats::na.omit(amr.ranges$dbeta)<=1) )
      stop("Malformed 'amr.ranges' object: 'dbeta' field is missing or is outside the range c(0,1)")
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
    .addEpimutations()
  }

  data.ranges <- GenomicRanges::granges(template.ranges)
  GenomicRanges::mcols(data.ranges) <- random.betas

  return(data.ranges)
}
