#' Search for aberrantly methylated regions
#'
#' @description
#' `getAMR` returns a `GRanges` object with all the aberrantly methylated
#' regions (AMRs) for all samples in a dataset.
#'
#' @details
#' In the provided dataset `getAMR` compares methylation beta values of each
#' sample with other samples to identify rare long-range methylation
#' aberrations.
#'
#' @param data.ranges A `GRanges` object with genomic locations and
#' corresponding beta values included as metadata.
#' @param data.samples A character vector with sample names (a subset of
#' metadata column names). If `NULL` (the default), then all samples (metadata
#' columns) are included in the analysis.
#' @param ramr.method A character scalar: when ramr.method is "IQR" (the
#' default), the filtering based on interquantile range is used (`iqr.cutoff`
#' value is then used as a threshold). When "beta" or "wbeta" - filtering based
#' on fitting non-weighted (EnvStats::ebeta) or weighted (ExtDist::eBeta) beta
#' distributions, respectively, is used, and `pval.cutoff` or `qval.cutoff` (if
#' not `NULL`) is used as a threshold. For "wbeta", weights directly correlate
#' with bin contents (number of values per bin) and inversly - with the
#' distances from the median value, thus narrowing the estimated distribution
#' and emphasizing outliers.
#' @param iqr.cutoff A single integer >= 1. Methylation beta values differing
#' from the median value by more than `iqr.cutoff` interquartile ranges are
#' considered to be significant (the default: 5).
#' @param pval.cutoff A numeric scalar (the default: 5e-2). Bonferroni
#' correction of `pval.cutoff` by the length of the `data.ranges` object is
#' used to calculate `qval.cutoff` if the latter is `NULL`.
#' @param qval.cutoff A numeric scalar. Used as a threshold for filtering based
#' on fitting non-weighted or weighted beta distributions: all p-values lower
#' than `qval.cutoff` are considered to be significant. If `NULL` (the default),
#' it is calculated using `pval.cutoff`
#' @param merge.window A positive integer. All significant (survived the
#' filtering stage) `data.ranges` genomic locations within this distance will be
#' merged to create AMRs (the default: 300).
#' @param min.cpgs A single integer >= 1. All AMRs containing less than
#' `min.cpgs` significant genomic locations are filtered out (the default: 7).
#' @param min.width A single integer >= 1 (the default). Only AMRs with the
#' width of at least `min.width` are returned.
#' @param exclude.range A numeric vector of length two. If not `NULL` (the
#' default), all `data.ranges` genomic locations with their median methylation
#' beta value within the `exclude.range` interval are filtered out.
#' @param cores A single integer >= 1. Number of processes for parallel
#' computation (the default: all but one cores).
#' @param ... Further arguments to be passed to `EnvStats::ebeta` or
#' `ExtDist::eBeta` functions.
#' @return The output is a `GRanges` object that contains all the aberrantly
#' methylated regions (AMRs) for all `data.samples` samples in `data.ranges`
#' object. The following metadata columns may be present:
#' \itemize{
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
#' \code{\link{simulateData}} for the generation of simulated test datasets,
#' and `ramr` vignettes for the description of usage and sample data.
#' @examples
#'   data(ramr)
#'   getAMR(ramr.data, ramr.samples, ramr.method="beta",
#'          min.cpgs=5, merge.window=1000, qval.cutoff=1e-3, cores=2)
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom GenomicRanges mcols reduce
#' @importFrom EnvStats ebeta
#' @importFrom ExtDist eBeta pBeta
#' @importFrom matrixStats rowMedians rowIQRs
#' @importFrom foreach foreach %dopar%
#' @importFrom methods as is
#' @importFrom stats median pbeta
#' @export
getAMR <- function (data.ranges,
                    data.samples=NULL,
                    ramr.method="IQR",
                    iqr.cutoff=5,
                    pval.cutoff=5e-2,
                    qval.cutoff=NULL,
                    merge.window=300,
                    min.cpgs=7,
                    min.width=1,
                    exclude.range=NULL, #c(0.3,0.7)
                    cores=max(1,parallel::detectCores()-1),
                    ...)
{
  if (!methods::is(data.ranges,"GRanges"))
    stop("'data.ranges' must be a GRanges object")
  if (is.null(data.samples))
    data.samples <- colnames(GenomicRanges::mcols(data.ranges))
  if (!all(data.samples %in% colnames(GenomicRanges::mcols(data.ranges))))
    stop("'data.ranges' metadata must include 'data.samples'")
  if (length(data.samples)<3)
    stop("at least three 'data.samples' must be provided")

  #####################################################################################
  getPValues.beta <- function (data.chunk, ...) {
    chunk.filt <- apply(data.chunk, 1, function (x) {
      x.median    <- stats::median(x, na.rm=TRUE)
      x[is.na(x)] <- x.median
      beta.fit <- suppressWarnings( EnvStats::ebeta(as.numeric(x), ...) )
      pvals    <- stats::pbeta(x, beta.fit$parameters[1], beta.fit$parameters[2])
      pvals[x>x.median] <- 1 - pvals[x>x.median]
      return(pvals)
    })
    return(t(chunk.filt))
  }
  getPValues.wbeta<- function (data.chunk, ...) {
    chunk.filt  <- apply(data.chunk, 1, function (x) {
      x.median    <- stats::median(x, na.rm=TRUE)
      x[is.na(x)] <- x.median
      # weight directly correlates with bin contents (number of values per bin)
      # and inversly - with the distance from the median value, thus narrowing
      # the estimated distribution and emphasizing outliers
      c           <- cut(x, c(0:100)/100)
      b           <- table(c)
      w           <- as.numeric(b[c]) * (1 - abs(x-x.median))
      beta.fit    <- suppressWarnings( ExtDist::eBeta(as.numeric(x), w, ...) )
      pvals       <- ExtDist::pBeta(x, params=beta.fit)
      pvals[x>x.median] <- 1 - pvals[x>x.median]
      return(pvals)
    })
    return(t(chunk.filt))
  }

  #####################################################################################


  doParallel::registerDoParallel(cores)
  cl <- parallel::makeCluster(cores)

  universe      <- getUniverse(data.ranges, merge.window=merge.window, min.cpgs=min.cpgs, min.width=min.width)
  universe.cpgs <- unlist(universe$revmap)

  betas <- as.matrix(mcols(data.ranges)[universe.cpgs,data.samples,drop=FALSE])
  if (is.null(qval.cutoff))
    qval.cutoff <- pval.cutoff/nrow(betas)

  chunks  <- split(seq_len(nrow(betas)), if (cores>1) cut(seq_len(nrow(betas)),cores) else 1)
  medians <- foreach (chunk=chunks, .combine=c) %dopar% matrixStats::rowMedians(betas[chunk,], na.rm=TRUE)

  if (ramr.method=="IQR") {
    iqrs <- foreach (chunk=chunks, .combine=c) %dopar% matrixStats::rowIQRs(betas[chunk,], na.rm=TRUE)
    betas.filtered <- (betas-medians)/iqrs
    betas.filtered[abs(betas.filtered)<iqr.cutoff]  <- NA
  } else if (ramr.method=="beta") {
    # multi-threaded EnvStats::ebeta (speed: mme=mmue>mle>>>fitdistrplus::fitdist)
    betas.filtered <- foreach (chunk=chunks) %dopar% getPValues.beta(betas[chunk,], ...)
    betas.filtered <- do.call(rbind, betas.filtered)
    betas.filtered[betas.filtered>=qval.cutoff] <- NA
  } else if (ramr.method=="wbeta") {
    betas.filtered <- foreach (chunk=chunks) %dopar% getPValues.wbeta(betas[chunk,], ...)
    betas.filtered <- do.call(rbind, betas.filtered)
    betas.filtered[betas.filtered>=qval.cutoff] <- NA
  } else {
    stop("unknown 'ramr.method'")
  }

  if (!is.null(exclude.range))
    betas.filtered[medians>=exclude.range[1] & medians<=exclude.range[2]] <- NA


  getMergedRanges <- function (column) {
    not.na <- which(!is.na(betas.filtered[,column]))
    ranges <- GenomicRanges::reduce(data.ranges[universe.cpgs[not.na]], min.gapwidth=merge.window, with.revmap=TRUE)
    if (length(ranges)>0) {
      ranges$ncpg   <- unlist(lapply(ranges$revmap, length))
      ranges$sample <- column
      ranges        <- subset(ranges, `ncpg`>=min.cpgs & `width`>=max(1,min.width))
      ranges$dbeta  <- vapply(ranges$revmap, function (revmap) {
        mean(betas[not.na[revmap],column,drop=FALSE] - medians[not.na[revmap]])
      }, numeric(1))

      if (ramr.method=="IQR") {
        ranges$xiqr   <- vapply(ranges$revmap, function (revmap) {
          mean(betas.filtered[not.na[revmap],column,drop=FALSE], na.rm=TRUE)
        }, numeric(1))
      }

      if (ramr.method=="beta" | ramr.method=="wbeta") {
        ranges$pval <- vapply(ranges$revmap, function (revmap) {
          return( 10**mean(log10(betas.filtered[not.na[revmap],column] + .Machine$double.xmin), na.rm=TRUE) )
        }, numeric(1))
      }

      ranges$revmap <- lapply(ranges$revmap, function (i) {universe.cpgs[not.na[i]]})
    }
    return(ranges)
  }

  amr.ranges <- foreach (column=colnames(betas.filtered)) %dopar% getMergedRanges(column)

  parallel::stopCluster(cl)
  return(unlist(methods::as(amr.ranges, "GRangesList")))
}
