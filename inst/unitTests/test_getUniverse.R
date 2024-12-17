test_getUniverse <- function () {
  data(ramr)
  RUnit::checkException(
    getUniverse(c())
  )
  RUnit::checkEquals(
    length( getUniverse(ramr.data, min.cpgs=1, merge.window=1) ),
    length(ramr.data)
  )
  RUnit::checkTrue(
    all( GenomicRanges::width( getUniverse(ramr.data, min.cpgs=10, merge.window=1000) ) > 1000 )
  )
}
