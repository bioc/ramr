test_getAMR <- function () {
  data(ramr)
  
  RUnit::checkException(
    getAMR(data.ranges=c())
  )
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, data.samples=ramr.samples[1:2])
  )
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, data.samples=c("a","b","c"))
  )
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, data.samples=ramr.samples, compute="zzz")
  )
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, compute.estimate="nmle")
  )
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, combine="comb-p")
  )
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, data.coverage=1)
  )
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, data.coverage=data.frame())
  )
  
  test.coverage <- as.data.frame(matrix(data=rbinom(n=3000*100, 100, 1:100/100), ncol=100))
  RUnit::checkException(
    getAMR(data.ranges=ramr.data, data.coverage=test.coverage[-1,])
  )
  
  amr.iqr.1 <- getAMR(data.ranges=ramr.data, compute="IQR", combine.min.cpgs=5, combine.window=10000, combine.threshold=5) #, ncores=1)
  amr.iqr.2 <- getAMR(data.ranges=ramr.data, data.coverage=test.coverage, compute="IQR", combine.min.cpgs=5, combine.window=10000, combine.threshold=5) #, ncores=2)
  RUnit::checkIdentical(
    amr.iqr.1,
    amr.iqr.2
  )
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.iqr.2, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.iqr.2, ramr.tp.nonunique))),
    c(6, 45)
  )
  
  amr.beta <- getAMR(data.ranges=ramr.data, data.samples=ramr.samples, data.coverage=test.coverage, compute="beta+binom", combine.min.cpgs=5, combine.window=10000, combine.threshold=1e-3)
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.beta, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.beta, ramr.tp.nonunique))),
    c(6, 45)
  )
  
  amr.wbeta <- getAMR(data.ranges=ramr.data, data.samples=ramr.samples, compute="beta+binom", compute.estimate="amle", compute.weights="sqrtInvDist",
                      combine.min.cpgs=5, combine.window=10000, combine.threshold=1e-5)
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.wbeta, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.wbeta, ramr.tp.nonunique))),
    c(6, 45)
  )
  
  amr.iqr.3 <- getAMR(data.ranges=ramr.data, compute="IQR", combine.min.cpgs=5, combine.window=10000, combine.threshold=5, exclude.range=c(0.1,0.9))
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.iqr.3, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.iqr.3, ramr.tp.nonunique))),
    c(2, 18)
  )
  
  
  ### tests to cover 0/1/NA values during fitting
  
  data.mcols <- GenomicRanges::mcols(ramr.data)[, -1]
  set.seed(1)
  data.mcols[sample(x=seq_len(nrow(data.mcols)), size=10), sample(x=seq_len(ncol(data.mcols)), size=10)] <- NA
  data.mcols[sample(x=seq_len(nrow(data.mcols)), size=10), sample(x=seq_len(ncol(data.mcols)), size=10)] <- 0
  data.mcols[sample(x=seq_len(nrow(data.mcols)), size=10), sample(x=seq_len(ncol(data.mcols)), size=10)] <- 1
  data.mcols[1500, ] <- NA
  mod.data <- GenomicRanges::granges(ramr.data)
  GenomicRanges::strand(mod.data) <- c("+", "-")
  GenomicRanges::mcols(mod.data) <- data.mcols
  mod.coverage <- as.data.frame(matrix(data=rbinom(n=nrow(data.mcols)*ncol(data.mcols), ncol(data.mcols), 1:100/100), ncol=ncol(data.mcols)))
  
  getAMR(data.ranges=mod.data, data.coverage=mod.coverage, compute="beta+binom", combine.min.cpgs=5, combine.window=10000, combine.threshold=1e-3, combine.ignore.strand=TRUE)
  getAMR(data.ranges=mod.data, compute="IQR", combine.min.cpgs=5, combine.window=10000, combine.threshold=3, combine.ignore.strand=TRUE)
  for (e in c("mom", "amle"))
    for (w in c("equal", "logInvDist", "sqrtInvDist", "invDist"))
      getAMR(data.ranges=mod.data, compute="beta+binom", compute.estimate=e, compute.weights=w, combine.min.cpgs=5, combine.window=10000, combine.threshold=1e-3)
  getAMR(data.ranges=mod.data, transform="linear", compute="beta+binom", compute.estimate=e, compute.weights=w, combine.min.cpgs=5, combine.window=10000, combine.threshold=1e-3)
}
