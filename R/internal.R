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
#' @importFrom utils head tail
## #' @importFrom Rcpp sourceCpp
## #' @useDynLib ramr, .registration=TRUE


# internal globals, constants and helper functions 
#

################################################################################
# Globals, unload
################################################################################

utils::globalVariables(c(
  "chunk", "column", "ncpg", "width", "..data.samples", ":=", "alpha", "color",
  "size", "start"
))

# .onUnload <- function (libpath) {library.dynam.unload("ramr", libpath)}

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
