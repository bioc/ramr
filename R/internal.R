## #' @importFrom data.table data.table
#' @importFrom Rcpp sourceCpp
#' @useDynLib ramr, .registration=TRUE


# internal globals, constants and helper functions 
#

################################################################################
# Globals, unload
################################################################################

utils::globalVariables(
  c("chunk", "column", "ncpg", "width")
)

.onUnload <- function (libpath) {library.dynam.unload("ramr", libpath)}

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
