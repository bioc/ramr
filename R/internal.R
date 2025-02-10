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

  fn <- paste("rcpp_prepare_data", transform, sep="_")
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

# descr: computes and filters by IQR
# value: void

.getAMR.IQR <- function (data.list,
                         threshold,
                         verbose)
{
  if (verbose) message("Computing IQR ", appendLF=FALSE)
  tm <- proc.time()

  rcpp_get_iqr(data=data.list)
  rcpp_compute_xiqr(data=data.list)
  rcpp_filter_threshold_xiqr(data=data.list, thr=threshold)

  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
}

################################################################################

# descr: fits beta distribution and filters by p-value
# value: void

.getAMR.beta <- function (data.list,
                          estimate,
                          weights,
                          coverage,
                          threshold,
                          verbose)
{
  if (verbose) message("Fitting beta distribution ", appendLF=FALSE)
  tm <- proc.time()

  fn.mean <- paste(
    "rcpp_get_meanvar",
    if (estimate=="mom") "ari" else "geo",
    weights,
    sep="_"
  )
  do.call(what=fn.mean, args=list(data=data.list))

  fn.fit <- paste("rcpp_fit_beta", estimate, sep="_")
  do.call(what=fn.fit, args=list(data=data.list))

  if (coverage==TRUE)
    rcpp_fit_binom(data=data.list)

  fn.logp <- paste0("rcpp_compute_logp_beta", if (coverage) "_binom")
  do.call(what=fn.logp, args=list(data=data.list))

  rcpp_filter_threshold_logp(data=data.list, thr=threshold)

  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
}

################################################################################

# descr: creates GenomicRanges
# value: GRanges object

.createGranges <- function (data.list,
                            compute,
                            window,
                            min.cpgs,
                            min.width,
                            ignore.strand,
                            verbose)
{
  if (verbose) message("Creating genomic ranges ", appendLF=FALSE)
  tm <- proc.time()

  fn.ranges <- paste(
    "rcpp_create_granges",
    if (ignore.strand) "unstranded" else "stranded",
    if (compute=="IQR") "xiqr" else "logp",
    sep="_"
  )
  amr.list <- do.call(what=fn.ranges, args=list(
    data=data.list, window=window,min_ncpg=min.cpgs, min_width=min.width
  ))

  amr.ranges <- GenomicRanges::GRanges(
    seqnames=amr.list$seqnames,
    ranges=IRanges::IRanges(start=amr.list$start, end=amr.list$end),
    strand=amr.list$strand,
    revmap=S4Vectors::I(amr.list$revmap),
    S4Vectors::DataFrame(
      amr.list[c("ncpg", "sample", "dbeta",
                 if (compute=="IQR") "xiqr" else "pval")]
    )
  )

  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(amr.ranges)
}

################################################################################

# descr: produces random data
# value: data.table

.getRandomValues <- function (data.list,
                              estimate,
                              weights,
                              nsamples,
                              sample.names,
                              verbose)
{
  if (verbose) message("Simulating data", appendLF=FALSE)
  tm <- proc.time()

  fn.mean <- paste(
    "rcpp_get_meanvar",
    if (estimate=="mom") "ari" else "geo",
    weights,
    sep="_"
  )
  do.call(what=fn.mean, args=list(data=data.list))

  fn.fit <- paste("rcpp_fit_beta", estimate, sep="_")
  do.call(what=fn.fit, args=list(data=data.list))

  random.values <- rcpp_generate_random_values(data=data.list)

  colnames(random.values) <- sample.names

  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(random.values)
}

################################################################################

# descr: introduces epimutations in random data
# value: data.table

.addEpimutations <- function (random.data,
                              amr.ranges,
                              verbose)
{
  if (verbose) message("Introducing epimutations", appendLF=FALSE)
  tm <- proc.time()

  amr.mcols <- data.frame(GenomicRanges::mcols(amr.ranges))
  for (i in seq_len(nrow(amr.mcols))) {
    revmap <- unlist(amr.mcols[i,"revmap"])
    dbeta  <- sign(0.5 - mean(random.betas[revmap,], na.omit=TRUE)) * amr.mcols[i,"dbeta"]
    sample <- amr.mcols[i,"sample"]
    random.betas[revmap, sample] <- random.betas[revmap, sample] + dbeta
  }

  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
}

################################################################################


