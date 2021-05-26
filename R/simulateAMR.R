#' @importFrom methods is
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
      v <- c(tail(v, -shift), head(v, shift))
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
  amr.ranges$dbeta  <- dbeta
  
  return(amr.ranges)
}