test_getAMR.obsolete <- function () {
  data(ramr)

  RUnit::checkException(
    getAMR.obsolete(c())
  )
  RUnit::checkException(
    getAMR.obsolete(ramr.data, ramr.samples[1:2])
  )
  RUnit::checkException(
    getAMR.obsolete(ramr.data, c("a","b","c"))
  )
  RUnit::checkException(
    getAMR.obsolete(ramr.data, ramr.samples, ramr.method="zzz")
  )

  amr.iqr.1 <- getAMR.obsolete(ramr.data, ramr.method="IQR", min.cpgs=5, merge.window=10000, iqr.cutoff=5, cores=1)
  amr.iqr.2 <- getAMR.obsolete(ramr.data, ramr.method="IQR", min.cpgs=5, merge.window=10000, iqr.cutoff=5, cores=2)
  RUnit::checkIdentical(
    amr.iqr.1,
    amr.iqr.2
  )
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.iqr.2, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.iqr.2, ramr.tp.nonunique))),
    c(6, 45)
  )

  amr.beta <- getAMR.obsolete(ramr.data, ramr.samples, ramr.method="beta", min.cpgs=5, merge.window=10000, qval.cutoff=1e-3, cores=2)
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.beta, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.beta, ramr.tp.nonunique))),
    c(6, 45)
  )

  amr.wbeta <- getAMR.obsolete(ramr.data, ramr.samples, ramr.method="wbeta", min.cpgs=5, merge.window=10000, qval.cutoff=1e-4, cores=2)
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.wbeta, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.wbeta, ramr.tp.nonunique))),
    c(6, 45)
  )

  amr.iqr.3 <- getAMR.obsolete(ramr.data, ramr.samples, ramr.method="IQR", min.cpgs=5, merge.window=10000, cores=2, exclude.range=c(0.1,0.9))
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.iqr.3, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.iqr.3, ramr.tp.nonunique))),
    c(2, 18)
  )

  amr.beinf <- getAMR.obsolete(ramr.data, ramr.samples, ramr.method="beinf", min.cpgs=5, merge.window=10000, qval.cutoff=1e-3, cores=2)
  RUnit::checkEquals(
    c(sum(GenomicRanges::countOverlaps(amr.beinf, ramr.tp.unique)), sum(GenomicRanges::countOverlaps(amr.beinf, ramr.tp.nonunique))),
    c(6, 45)
  )

  ### if suggested library is not available
  test.env <- new.env()
  assign(x="is.test.environment", value=TRUE, envir=test.env)
  test.func <- function(f, env, ...) {
    environment(f) <- env
    f(...)
  }
  RUnit::checkException(
    test.func(f=getAMR.obsolete, env=test.env,
              ramr.data, ramr.samples, ramr.method="beta", min.cpgs=5, merge.window=10000, qval.cutoff=1e-3, cores=2)
  )
}
