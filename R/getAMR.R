getAMR <- function (data.ranges,
                    data.samples,
                    ramr.method="IQR",
                    iqr.cutoff=5,
                    pval.cutoff=5e-2,
                    qval.cutoff=NULL,
                    merge.window=300,
                    min.cpgs=7,
                    min.width=NULL,
                    exclude.range=NULL, #c(0.3,0.7)
                    cores=max(1,detectCores()-1),
                    ...)
{
  if (class(data.ranges)!="GRanges")
    stop("'data.ranges' must be be a GRanges object")
  if (!all(data.samples %in% colnames(mcols(data.ranges))))
    stop("'data.ranges' metadata must include 'data.samples'")
  if (length(data.samples)<3)
    stop("at least three 'data.samples' must be provided")

  #####################################################################################
  getPValues.beta <- function(data.chunk, ...) {
    chunk.filt <- apply(data.chunk, 1, function (x) {
      x.median    <- median(x, na.rm=TRUE)
      x[is.na(x)] <- x.median
      # x.pvals <- tryCatch(
      #   {
      #     beta.fit <- suppressWarnings( EnvStats::ebeta(as.numeric(x), ...) )
      #     pvals    <- pbeta(x, beta.fit$parameters[1], beta.fit$parameters[2])
      #     pvals[x>x.median] <- 1 - pvals[x>x.median]
      #     pvals
      #   },
      #   error   = function (e) {
      #     return(rep(1, length(x)))
      #   }
      # )
      # return(x.pvals)
      beta.fit <- suppressWarnings( EnvStats::ebeta(as.numeric(x), ...) )
      pvals    <- pbeta(x, beta.fit$parameters[1], beta.fit$parameters[2])
      pvals[x>x.median] <- 1 - pvals[x>x.median]
      return(pvals)
    })
    return(t(chunk.filt))
  }
  getPValues.wbeta<- function(data.chunk, ...){
    chunk.filt  <- apply(data.chunk, 1, function(x){
      x.median    <- median(x, na.rm=TRUE)
      x[is.na(x)] <- x.median
      # # weight: bin contents
      #   c           <- cut(x, c(0:100)/100)
      #   b           <- table(c)
      #   w           <- as.numeric(b[c])
      # # weight: distance from median
      #   w           <- (1 - abs(x-x.median))**3
      # weight: bin contents and distance
      c           <- cut(x, c(0:100)/100)
      b           <- table(c)
      w           <- as.numeric(b[c]) * (1 - abs(x-x.median))
      beta.fit    <- suppressWarnings( ExtDist::eBeta(as.numeric(x), w, ...) )
      pvals       <- pBeta(x, params=beta.fit)
      pvals[x>x.median] <- 1 - pvals[x>x.median]
      return(pvals)
    })
    return(t(chunk.filt))
  }
  getPValues.logn <- function(data.chunk, ...){
    chunk.filt  <- apply(data.chunk, 1, function(x){
      x.quantiles <- quantile(x, na.rm=TRUE)
      x[is.na(x)] <- x.quantiles[3]
      logn.fit    <- suppressWarnings( logitnorm::twCoefLogitnorm(x.quantiles[3], x.quantiles[4], perc=0.75, ...) )
      pvals       <- plogitnorm(x, logn.fit[1], logn.fit[2])
      pvals[x>x.quantiles[3]] <- 1 - pvals[x>x.quantiles[3]]
      return(pvals)
    })
    return(t(chunk.filt))
  }
  #####################################################################################


  registerDoParallel(cores)
  cl <- makeCluster(cores)

  # TODO: all multicore                 - DONE (almost)
  # TODO: tile window                   -
  # TODO: distributions for outliers    - DONE (beta & logitnorm, but maybe try other skew-normal such as fGarch::snormFit)

  universe      <- getUniverse(data.ranges, merge.window=merge.window, min.cpgs=min.cpgs)
  universe.cpgs <- unlist(universe$revmap)

  betas   <- as.matrix(mcols(data.ranges)[universe.cpgs,data.samples,drop=FALSE])
  if (is.null(qval.cutoff))
    qval.cutoff <- pval.cutoff/nrow(betas)
  # names(dimnames(betas)) <- c("cpg", "sample")

  chunks <- split(1:nrow(betas), cut(1:nrow(betas),cores))
  medians <- foreach (chunk=chunks, .combine=c) %dopar% matrixStats::rowMedians(betas[chunk,], na.rm=TRUE)

  if (!is.null(exclude.range))
    medians[medians %between% exclude.range] <- NA

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
  } else if (ramr.method=="logn") {
    warning("logit-normal ramr.method doesn't work really...")
    # multi-threaded logitnorm::twCoefLogitnorm
    betas.filtered <- foreach (chunk=chunks) %dopar% getPValues.logn(betas[chunk,], ...)
    betas.filtered <- do.call(rbind, betas.filtered)
    betas.filtered[betas.filtered>=qval.cutoff] <- NA
  } else {
    stop("unknown 'ramr.method'")
  }

  getMergedRanges <- function (column) {
    not.na        <- which(!is.na(betas.filtered[,column]))
    ranges        <- GenomicRanges::reduce(data.ranges[universe.cpgs[not.na]], min.gapwidth=merge.window, with.revmap=TRUE)
    if (length(ranges)>0) {
      ranges$ncpg   <- unlist(lapply(ranges$revmap, length))
      ranges$sample <- column
      ranges        <- subset(ranges, ncpg>=min.cpgs & width>=max(3,min.width))
      ranges$dbeta  <- sapply(ranges$revmap, function (revmap) {
        mean(betas[not.na[revmap],column,drop=FALSE] - medians[not.na[revmap]])
      })

      if (ramr.method=="IQR") {
        ranges$xiqr   <- sapply(ranges$revmap, function (revmap) {
          mean(betas.filtered[not.na[revmap],column,drop=FALSE], na.rm=TRUE)
        })
        # } else {
        #   ranges$xiqr   <- sapply(ranges$revmap, function (revmap) {
        #     xiqrs <- (betas[not.na[revmap],column,drop=FALSE] - medians[not.na[revmap]]) / matrixStats::rowIQRs(betas[not.na[revmap],,drop=FALSE], na.rm=TRUE)
        #     mean(xiqrs)
        #   })
      }

      if (ramr.method=="beta") {
        ranges$pval <- sapply(ranges$revmap, function (revmap) {
          return( 10**mean(log10(betas.filtered[not.na[revmap],column] + .Machine$double.xmin), na.rm=TRUE) )
        })
        # } else {
        #   ranges$pval.beta <- sapply(ranges$revmap, function (revmap) {
        #     range.pvals <- getPValues.beta(betas[not.na[revmap],,drop=FALSE], ...)
        #     range.pvals[range.pvals==0] <- .Machine$double.xmin
        #     return( 10**mean(log10(range.pvals[,column]), na.rm=TRUE) )
        #   })
      }

      if (ramr.method=="wbeta") {
        ranges$pval <- sapply(ranges$revmap, function (revmap) {
          return( 10**mean(log10(betas.filtered[not.na[revmap],column] + .Machine$double.xmin), na.rm=TRUE) )
        })
        # } else {
        #   ranges$pval.wbeta <- sapply(ranges$revmap, function (revmap) {
        #     range.pvals <- getPValues.beta(betas[not.na[revmap],,drop=FALSE], ...)
        #     range.pvals[range.pvals==0] <- .Machine$double.xmin
        #     return( 10**mean(log10(range.pvals[,column]), na.rm=TRUE) )
        #   })
      }

      if (ramr.method=="logn") {
        ranges$pval <- sapply(ranges$revmap, function (revmap) {
          return( 10**mean(log10(betas.filtered[not.na[revmap],column] + .Machine$double.xmin), na.rm=TRUE) )
        })
        # } else {
        #   ranges$pval.logn <- sapply(ranges$revmap, function (revmap) {
        #     range.pvals <- getPValues.logn(betas[not.na[revmap],,drop=FALSE], ...)
        #     range.pvals[range.pvals==0] <- .Machine$double.xmin
        #     return( 10**mean(log10(range.pvals[,column]), na.rm=TRUE) )
        #   })
      }

      ranges$revmap <- lapply(ranges$revmap, function (i) {universe.cpgs[not.na[i]]})
    }
    return(ranges)
  }

  amr.ranges <- foreach (column=colnames(betas.filtered)) %dopar% getMergedRanges(column)

  stopCluster(cl)
  return(unlist(as(amr.ranges, "GRangesList")))
}
