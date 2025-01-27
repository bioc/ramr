#' @importFrom BiocGenerics relist
#' @importFrom data.table as.data.table melt.data.table
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng%
#' @importFrom EnvStats ebeta
#' @importFrom ExtDist eBeta pBeta
#' @importFrom foreach foreach
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist pBEINF
#' @importFrom GenomicRanges mcols `mcols<-` granges reduce findOverlaps
#' @importFrom IRanges subsetByOverlaps
#' @importFrom matrixStats rowMedians rowIQRs
#' @importFrom methods as is
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom S4Vectors queryHits
#' @importFrom stats median na.omit rbeta pbeta setNames
#' @importFrom utils head tail packageVersion
#' @importFrom Rcpp sourceCpp
#' @useDynLib ramr, .registration=TRUE


# internal globals, constants and helper functions
#

################################################################################
# Globals, unload, attach
################################################################################

utils::globalVariables(c(
  "chunk", "column", "ncpg", "width", "..data.samples", ":=", "alpha", "color",
  "size", "start"
))

.onUnload <- function (libpath) {library.dynam.unload("ramr", libpath)}

.onAttach <- function(libname, pkgname) {
  if(interactive()) {
    max.threads <- rcpp_test_omp()
    msg <- ifelse(
      max.threads<0,
      paste0(
        "Multithreading (OpenMP) is not available.\n",
        "Check how to enable it at https://github.com/BBCG/ramr"
      ),
      sprintf(
        "ramr v%s using %i out of %i available threads",
        utils::packageVersion("ramr"), max(1, max.threads %/% 2), max.threads
      )
    )
    packageStartupMessage(msg)
  }
  invisible()
}

################################################################################
# Constants
################################################################################

# descr: ...

#

################################################################################
# Functions
################################################################################

# descr: preprocesses data
# value: list

.preprocessData <- function (data.ranges,
                             data.samples,
                             data.coverage,
                             transform,
                             exclude.range,
                             ncores,
                             verbose)
{
  if (verbose) message("Preprocessing data ", appendLF=FALSE)
  tm <- proc.time()
  
  fn <- paste0("rcpp_prepare_data_", transform)
  data.object <- do.call(what=fn, args=list(
    seqnames=S4Vectors::runValue(GenomeInfoDb::seqnames(data.ranges)),
    seqrunlens=S4Vectors::runLength(GenomeInfoDb::seqnames(data.ranges)),
    start=BiocGenerics::start(data.ranges),
    strand=as.factor(BiocGenerics::strand(data.ranges)),
    mcols=as.data.frame(GenomicRanges::mcols(data.ranges), optional=TRUE),
    coverage=data.coverage,
    exclude_lower=exclude.range[1],
    exclude_upper=exclude.range[2]
  ))
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(data.object)
}

################################################################################


