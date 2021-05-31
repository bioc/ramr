#' Simulate a set of aberrantly methylated regions
#'
#' @description
#' `simulateAMR` returns a `GRanges` object containing a set of randomly
#' selected aberrantly methylated regions (AMRs) to be used as an input for the
#' `simulateData` method.
#'
#' @details
#' Using provided template (`GRanges` object) `simulateAMR` randomly selects
#' genomic regions satisfying various criteria (number of CpGs, width of the
#' region) and assigns them to samples according to specified parameters
#' (number of AMRs per sample, number of samples per AMR). Its output is meant
#' to be used as the set of true positive AMRs for the `simulateData` function.
#'
#' @param template.ranges A `GRanges` object with genomic locations (same object
#' must be supplied to this and to the `simulateData` functions).
#' @param nsamples A single integer >= 1 indicating the number of samples to
#' which AMRs will be assigned.
#' @param exclude.ranges A `GRanges` object with genomic locations. None of the
#' simulated AMRs in the output will overlap with any of regions from
#' `exclude.ranges`. If `NULL` (the default), AMRs are not restricted by their
#' genomic location.
#' @param regions.per.sample A single integer >= 1 (the default). Number of AMRs
#' to be assigned to every sample. Message is shown and the `regions.per.sample`
#' value is limited to `maxnAMR %/% nsamples` if it is greater than this number
#' (where `maxnAMR` is the maximum number of potential AMRs for the
#' `template.ranges`).
#' @param samples.per.region A single integer >= 1 (the default). Number of
#' samples to which the same AMR will be assigned. Message is shown and the
#' `samples.per.region` value is limited to `nsamples` if the former is greater
#' than the latter.
#' @param sample.names A character vector with sample names. If `NULL` (the
#' default), sample names will be computed as
#' `paste0("sample", seq_len(nsamples))`. When specified, the length of the
#' `sample.names` vector must not be smaller than the value of `nsamples`.
#' @param merge.window A positive integer. All `template.ranges` genomic
#' locations within this distance will be merged to create a list of potential
#' AMRs (which will be later filtered from regions overlapping with any regions
#' from the `exclude.ranges`).
#' @param min.cpgs A single integer >= 1. All AMRs containing less than
#' `min.cpgs` genomic locations are filtered out. The default: 7.
#' @param max.cpgs A single integer >= 1. All AMRs containing more than
#' `max.cpgs` genomic locations are filtered out. The default: `Inf`.
#' @param min.width A single integer >= 1 (the default). Only AMRs with the
#' width of at least `min.width` are returned.
#' @param dbeta A single non-negative numeric value in the range [0,1] or a
#' numeric vector of such values (with as many elements as there are AMRs).
#' Used to populate the `dbeta` metadata column, defines a desired absolute
#' deviation of corresponding AMR from the median for the `simulateData`
#' function.
#' @return The output is a `GRanges` object that contains a subset of
#' aberrantly methylated regions (AMRs) randomly selected from all the possible
#' AMRs for the provided `template.ranges` object. The following metadata
#' columns are included:
#' \itemize{
#'   \item `revmap` -- integer list of `template.ranges` genomic locations that
#'   are included in this AMR region
#'   \item `ncpg` -- number of `template.ranges` genomic locations within
#'   this AMR region
#'   \item `sample` -- an identifier of a sample to which
#'   corresponding AMR belongs
#'   \item `dbeta` -- equals to supplied `dbeta` parameter
#' }
#' @seealso \code{\link{simulateData}} for the generation of simulated test
#' data sets, \code{\link{getAMR}} for identification of AMRs,
#' \code{\link{plotAMR}} for plotting AMRs, \code{\link{getUniverse}}
#' for info on enrichment analysis, and `ramr` vignettes for the description of
#' usage and sample data.
#' @examples
#'   data(ramr)
#'   amrs.unique <-
#'     simulateAMR(ramr.data, nsamples=4, regions.per.sample=2,
#'                 min.cpgs=5, merge.window=1000, dbeta=0.2)
#'   amrs.nonunique <-
#'     simulateAMR(ramr.data, nsamples=3, exclude.ranges=amrs.unique,
#'                 samples.per.region=2, min.cpgs=5, merge.window=1000)
#' @importFrom methods is
#' @importFrom utils head tail
#' @importFrom IRanges `%outside%`
#' @export
simulateAMR <- function (template.ranges,
                         nsamples,
                         exclude.ranges=NULL,
                         regions.per.sample=1,
                         samples.per.region=1,
                         sample.names=NULL,
                         merge.window=300,
                         min.cpgs=7,
                         max.cpgs=Inf,
                         min.width=1,
                         dbeta=0.25)
{
  if (!methods::is(template.ranges,"GRanges"))
    stop("'template.ranges' must be a GRanges object")
  if (!is.null(sample.names) & length(sample.names)<nsamples)
    stop("'sample.names' length must be greater or equal to 'nsamples'")

  #####################################################################################

  repRotate <- function(v, times=1, shift=1) {
    r <- v
    for (i in seq_len(times-1)) {
      v <- c(utils::tail(v, -shift), utils::head(v, shift))
      r <- c(r, v)
    }
    return(r)
  }

  #####################################################################################

  universe.ranges <- getUniverse(template.ranges, merge.window=merge.window, min.cpgs=min.cpgs, min.width=min.width)
  universe.ranges <- subset(universe.ranges, `ncpg`<=max.cpgs)
  if (methods::is(exclude.ranges,"GRanges"))
    universe.ranges <- subset(universe.ranges, universe.ranges %outside% exclude.ranges)
  if (is.null(sample.names))
    sample.names <- paste0("sample", seq_len(nsamples))

  if (nsamples > length(universe.ranges)) {
    nsamples <- length(universe.ranges)
    message("'nsamples' has been set to the maximum possible value of ", nsamples)
  }
  if (regions.per.sample > length(universe.ranges) %/% nsamples) {
    regions.per.sample <- length(universe.ranges) %/% nsamples
    message("'regions.per.sample' has been set to the maximum possible value of ", regions.per.sample)
  }
  if (samples.per.region > nsamples) {
    samples.per.region <- nsamples
    message("'samples.per.region' has been set to the maximum possible value of ", nsamples)
  }

  amr.ranges        <- rep(sample(universe.ranges, nsamples * regions.per.sample), samples.per.region)
  amr.ranges$sample <- repRotate(rep(sample(sample.names, nsamples), regions.per.sample), times=samples.per.region, shift=1)
  amr.ranges$dbeta  <- as.numeric(dbeta)

  return(amr.ranges)
}
