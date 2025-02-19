#' Search for aberrantly methylated regions
#'
#' @description
#' `getAMR` returns a `GRanges` object with aberrantly methylated
#' regions (AMRs / epimutations) for all samples in a data set.
#'
#' @details
#' In the provided data set, `getAMR` finds stretches of outlier beta values
#' to identify rare long-range methylation
#' aberrations (epimutations) in one or several samples. Other methods for
#' differential methylation analysis rely on between-group comparisons ---
#' `getAMR` performs this comparison within-sample, which is not only faster,
#' but also more sensitive. The logic of computations is described below.
#'
#' \subsection{Compute}{
#' Currently, only two methods are supported: ... ...
#' \subsection{IQR}{
#' here's IQR
#' }
#'
#' or
#'
#' \describe{
#'   \item{IQR}{
#'     When `compute=="IQR"`, for every genomic location (CpG) in
#'     `data.ranges` the IQR-normalized deviation from the median value is
#'     calculated, and all CpGs with such normalized deviation not smaller
#'     than the `iqr.cutoff` are retained.
#'   }
#' }
#'
#'  For
#' `ramr.method %in% c("beta", "wbeta", "beinf")`: parameters of beta
#' distribution are estimated by means of `EnvStats::ebeta` (beta distribution),
#' `ExtDist::eBeta` (weighted beta destribution), or `gamlss.dist::BEINF` (zero
#' and one inflated beta distribution) functions, respectively. These
#' parameters are then used to calculate the probability values, followed by the
#' filtering when all CpGs with p-values not greater than `qval.cutoff` are
#' retained. Another filtering is then performed to exclude all CpGs within
#' `exclude.range`. Next, the retained (significant) CpGs are merged within
#' the window of `merge.window`, and final filtering is applied to AMR genomic
#' ranges (by `min.cpgs` and `min.width`).
#' }
#'
#' \subsection{Combine}{
#' this and that
#' }
#'
#' @param data.ranges A `GRanges` object with genomic locations and
#' corresponding beta values included as metadata.
#' @param data.samples A character vector with sample names (a subset of
#' metadata column names). If `NULL` (the default), then all samples (metadata
#' columns) are included in the analysis.
#' @param data.coverage description
#' @param transform description
#' @param exclude.range A numeric vector of length two. If not `NULL` (the
#' default), all `data.ranges` genomic locations with their median methylation
#' beta value within the `exclude.range` interval are filtered out.
#' @param compute A character scalar: when ramr.method is "IQR" (the
#' default), the filtering based on interquantile range is used (`iqr.cutoff`
#' value is then used as a threshold). When "beta", "wbeta" or "beinf" -
#' filtering based on fitting non-weighted (`EnvStats::ebeta`), weighted
#' (`ExtDist::eBeta`) or zero-and-one inflated (`gamlss.dist::BEINF`)
#' beta distributions, respectively, is used, and `pval.cutoff` or `qval.cutoff`
#' (if not `NULL`) is used as a threshold. For "wbeta", weights directly
#' correlate with bin contents (number of values per bin) and inversly - with
#' the distances from the median value, thus narrowing the estimated
#' distribution and emphasizing outliers.
#' @param compute.estimate description
#' @param compute.weights description
#' @param combine description
#' @param combine.threshold A single integer >= 1. Methylation beta values differing
#' from the median value by more than `iqr.cutoff` interquartile ranges are
#' considered to be significant (the default: 5).
#' param pval.cutoff A numeric scalar (the default: 5e-2). Bonferroni
#' correction of `pval.cutoff` by the length of the `data.samples` object is
#' used to calculate `qval.cutoff` if the latter is `NULL`.
#' param qval.cutoff A numeric scalar. Used as a threshold for filtering based
#' on fitting non-weighted or weighted beta distributions: all p-values lower
#' than `qval.cutoff` are considered to be significant. If `NULL` (the default),
#' it is calculated using `pval.cutoff`
#' @param combine.window A positive integer. All significant (survived the
#' filtering stage) `data.ranges` genomic locations within this distance will be
#' merged to create AMRs (the default: 300).
#' @param combine.min.cpgs A single integer >= 1. All AMRs containing less than
#' `min.cpgs` significant genomic locations are filtered out (the default: 7).
#' @param combine.min.width A single integer >= 1 (the default). Only AMRs with the
#' width of at least `min.width` are returned.
#' @param combine.ignore.strand description
#' @param ncores A single integer >= 1. Number of processes for parallel
#' computation (the default: all but one cores). Results of parallel processing
#' are fully reproducible when the same seed is used (thanks to doRNG).
#' @param verbose Boolean to report progress and timings (default: TRUE).
#' @return The output is a `GRanges` object that contains all the aberrantly
#' methylated regions (AMRs) for all `data.samples` samples in `data.ranges`
#' object. The following metadata columns may be present:
#' \itemize{
#'   \item `revmap` -- integer list of significant CpGs (`data.ranges` genomic
#'   locations) that are included in this AMR region
#'   \item `ncpg` -- number of significant CpGs within this AMR region
#'   \item `sample` -- contains an identifier of a sample to which
#'   corresponding AMR belongs
#'   \item `dbeta` -- average deviation of beta values for significant CpGs from
#'   their corresponding median values
#'   \item `pval` -- geometric mean of p-values for significant CpGs
#'   \item `xiqr` -- average IQR-normalised deviation of beta values for
#'   significant CpGs from their corresponding median values
#' }
#' @seealso \code{\link{plotAMR}} for plotting AMRs, \code{\link{getUniverse}}
#' for info on enrichment analysis, \code{\link{simulateAMR}} and
#' \code{\link{simulateData}} for the generation of simulated test data sets,
#' and `ramr` vignettes for the description of usage and sample data.
#' @examples
#'   data(ramr)
#'   getAMR(data.ranges=ramr.data, data.samples=ramr.samples,
#'          compute="beta+binom", compute.estimate="amle",
#'          combine.min.cpgs=5, combine.window=1000, combine.threshold=1e-3)
#' @export
getAMR <- function (data.ranges,
                    data.samples=NULL,
                    data.coverage=NULL,
                    transform=c("identity", "linear"),
                    exclude.range=NULL,
                    compute=c("IQR", "beta+binom"),
                    compute.estimate=c("mom", "amle", "nmle"),
                    compute.weights=c("equal", "logInvDist", "sqrtInvDist", "invDist"),
                    combine=c("threshold", "comb-p"),
                    combine.threshold=ifelse(compute=="IQR", 5, 1e-3),
                    combine.window=300,
                    combine.min.cpgs=7,
                    combine.min.width=1,
                    combine.ignore.strand=FALSE,
                    ncores=NULL,
                    verbose=TRUE)
{
  if (!methods::is(data.ranges,"GRanges"))
    stop("'data.ranges' must be a GRanges object")
  data.mcols <- GenomicRanges::mcols(data.ranges)
  if (is.null(data.samples))
    data.samples <- colnames(data.mcols)
  if (!all(data.samples %in% colnames(data.mcols)))
    stop("'data.ranges' metadata must include 'data.samples'")
  if (!is.null(data.coverage) &
      (!methods::is(data.coverage,"data.frame") |
       !identical(dim(data.mcols), dim(data.coverage))))
    stop("When provided, 'data.coverage' must be a 'data.frame' object",
         " of the same dimensions as 'data.ranges' metadata")
  if (length(data.samples)<3)
    stop("at least three 'data.samples' must be provided")

  if (is.null(data.coverage))
    data.coverage <- data.frame()
  if (is.null(exclude.range))
    exclude.range <- c(2,0) # // <= than 2 and >= than 0
  transform <- match.arg(transform)
  compute <- match.arg(compute)
  compute.estimate <- match.arg(compute.estimate)
  compute.weights <- match.arg(compute.weights)
  combine <- match.arg(combine)

  if (compute.estimate=="nmle")
    stop("compute.estimate=='nmle' is not available yet")
  if (combine=="comb-p")
    stop("combine=='comb-p' is not available yet")

  #####################################################################################

  .data <- .preprocessData(
    data.ranges=data.ranges, data.samples=data.samples,
    data.coverage=data.coverage, transform=transform,
    exclude.range=exclude.range, ncores=ncores, verbose=verbose
  )

  if (compute=="IQR") {
    .getAMR.IQR(
      data.list=.data, threshold=combine.threshold, verbose=verbose
    )
  } else if (compute=="beta+binom") {
    .getAMR.beta(
      data.list=.data,
      estimate=compute.estimate,
      weights=compute.weights,
      coverage=identical(dim(data.mcols), dim(data.coverage)),
      threshold=log(combine.threshold),
      verbose=verbose
    )
  }

  amr.ranges <- .createGranges(
    data.list=.data,
    compute=compute,
    window=combine.window,
    min.cpgs=combine.min.cpgs,
    min.width=combine.min.width,
    ignore.strand=combine.ignore.strand,
    verbose=verbose
  )

  return(amr.ranges)
}
