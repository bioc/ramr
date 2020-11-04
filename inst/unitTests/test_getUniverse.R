test_getUniverse <- function () {
  checkEquals(
    length( getUniverse(ramr.data, min.cpgs=1, merge.window=1) ),
    length(ramr.data)
  )
  checkTrue(
    all( width( getUniverse(ramr.data, min.cpgs=10, merge.window=1000) ) > 1000 )
  )
}
