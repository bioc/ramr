test_getAMR <- function () {
  data(ramr)
  RUnit::checkEquals(
    length(plotAMR(ramr.data, ramr.samples, ramr.tp.nonunique)),
    length(reduce(ramr.tp.nonunique))
  )
}
