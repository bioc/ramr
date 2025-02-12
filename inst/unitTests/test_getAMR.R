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
  
  amr.iqr.1 <- getAMR(data.ranges=ramr.data, compute="IQR", combine.min.cpgs=5, combine.window=10000, combine.threshold=5) #, ncores=1)
  amr.iqr.2 <- getAMR(data.ranges=ramr.data, compute="IQR", combine.min.cpgs=5, combine.window=10000, combine.threshold=5) #, ncores=2)
  RUnit::checkIdentical(
    amr.iqr.1,
    amr.iqr.2
  )
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.iqr.2, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.iqr.2, ramr.tp.nonunique))),
    c(6, 45)
  )
  
  amr.beta <- getAMR(data.ranges=ramr.data, data.samples=ramr.samples, compute="beta+binom", combine.min.cpgs=5, combine.window=10000, combine.threshold=1e-3)
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
}
