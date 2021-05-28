#' Merges, filters and outputs all genomic regions of a given `GRanges` object
#'
#' @description
#' `getUniverse` returns a GRanges object with all the genomic regions in a
#' dataset, that can be used for AMR enrichment analysis
#'
#' @details
#' In the provided dataset `getUniverse` merges and outputs all the genomic
#' regions that satisfy filtering criteria, thus creating a `GRanges` object to
#' be used as a reference set of genomic regions for AMR enrichment analysis.
#'
#' @param data.ranges A `GRanges` object with genomic locations and
#' corresponding beta values included as metadata.
#' @param merge.window A positive integer. All `data.ranges` genomic locations
#' within this distance will be merged (the default: 300).
#' @param min.cpgs A single integer >= 1. All genomic regions containing less
#' than `min.cpgs` genomic locations are filtered out (the default: 7).
#' @param min.width A single integer >= 1 (the default). Only regions with the
#' width of at least `min.width` are returned.
#' @return The output is a `GRanges` object that contain all the genomic regions
#' in `data.ranges` object (in other words, all potential AMRs).
#' @seealso \code{\link{getAMR}} for identification of AMRs,
#' \code{\link{plotAMR}} for plotting AMRs, \code{\link{simulateAMR}} and
#' \code{\link{simulateData}} for the generation of simulated test datasets,
#' and `ramr` vignettes for the description of usage and sample data.
#' @examples
#'   data(ramr)
#'   universe <- getUniverse(ramr.data, min.cpgs=5, merge.window=1000)
#' \donttest{
#'   # identify AMRs
#'   amrs <- getAMR(ramr.data, ramr.samples, ramr.method="beta", min.cpgs=5,
#'                  merge.window=1000, qval.cutoff=1e-3, cores=2)
#'
#'   # AMR enrichment analysis using LOLA
#'   library(LOLA)
#'   # download LOLA region databases from http://databio.org/regiondb
#'   hg19.extdb.file <- system.file("LOLAExt", "hg19", package="LOLA")
#'   if (file.exists(hg19.extdb.file)) {
#'     hg19.extdb  <- loadRegionDB(hg19.extdb.file)
#'     runLOLA(amrs, universe, hg19.extdb, cores=1, redefineUserSets=TRUE)
#'   }
#' }
#' @importFrom GenomicRanges reduce
#' @importFrom methods is
#' @export
getUniverse <- function (data.ranges,
                         merge.window=300,
                         min.cpgs=7,
                         min.width=1)
{
  if (!methods::is(data.ranges,"GRanges"))
    stop("'data.ranges' must be a GRanges object")

  universe.ranges <- GenomicRanges::reduce(data.ranges, min.gapwidth=merge.window, with.revmap=TRUE)
  if (length(universe.ranges)>0) {
    universe.ranges$ncpg <- unlist(lapply(universe.ranges$revmap, length))
    universe.ranges      <- subset(universe.ranges, `ncpg`>=min.cpgs & `width`>=max(1,min.width))
  }
  return(universe.ranges)
}
