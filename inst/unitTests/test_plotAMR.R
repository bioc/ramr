test_plotAMR <- function () {
  data(ramr)
  RUnit::checkEquals(
    length(plotAMR(ramr.data, ramr.samples, ramr.tp.nonunique)),
    length(reduce(ramr.tp.nonunique))
  )
  RUnit::checkEquals(
    length(plotAMR(ramr.data, NULL, ramr.tp.unique)),
    length(reduce(ramr.tp.unique))
  )
}
