test_getAMR <- function () {
  data(ramr)
  RUnit::checkException(
    getAMR(c())
  )
  RUnit::checkException(
    getAMR(ramr.data, ramr.samples[1:2])
  )
  RUnit::checkException(
    getAMR(ramr.data, c("a","b","c"))
  )
  RUnit::checkException(
    getAMR(ramr.data, ramr.samples, ramr.method="zzz")
  )
  RUnit::checkEquals(
    length( getAMR(ramr.data, ramr.method="IQR", min.cpgs=5, merge.window=10000, iqr.cutoff=5, cores=1) ),
    length( c(ramr.tp.unique,ramr.tp.nonunique) )
  )
  RUnit::checkEquals(
    length( getAMR(ramr.data, ramr.samples, ramr.method="beta", min.cpgs=5, merge.window=10000, qval.cutoff=1e-3, cores=2) ),
    length( c(ramr.tp.unique,ramr.tp.nonunique) )
  )
  RUnit::checkEquals(
    length( getAMR(ramr.data, ramr.samples, ramr.method="wbeta", min.cpgs=5, merge.window=10000, qval.cutoff=1e-4, cores=2) ),
    length( c(ramr.tp.unique,ramr.tp.nonunique) )
  )
  RUnit::checkEquals(
    length( getAMR(ramr.data, ramr.samples, ramr.method="IQR", min.cpgs=5, merge.window=10000, qval.cutoff=1e-3, cores=2, exclude.range=c(0.1,0.9)) ),
    8
  )
}
