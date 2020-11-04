test_getAMR <- function () {
  checkEquals(
    length( getAMR(ramr.data, ramr.samples, ramr.method="IQR", min.cpgs=5, merge.window=10000, iqr.cutoff=5, cores=1) ),
    length( c(ramr.tp.unique,ramr.tp.nonunique) )
  )
  checkEquals(
    length( getAMR(ramr.data, ramr.samples, ramr.method="beta", min.cpgs=5, merge.window=10000, qval.cutoff=1e-3, cores=2) ),
    length( c(ramr.tp.unique,ramr.tp.nonunique) )
  )
  checkEquals(
    length( getAMR(ramr.data, ramr.samples, ramr.method="wbeta", min.cpgs=5, merge.window=10000, qval.cutoff=1e-4, cores=2) ),
    length( c(ramr.tp.unique,ramr.tp.nonunique) )
  )
}
