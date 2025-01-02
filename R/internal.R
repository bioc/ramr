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
#' @importFrom stats median na.omit rbeta pbeta
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
# Functions: ...
################################################################################

# descr: ...
# value: ...
