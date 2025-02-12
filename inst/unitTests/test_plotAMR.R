test_plotAMR <- function () {
  data(ramr)
  
  RUnit::checkException(
    plotAMR(data.ranges=sort(ramr.data, decreasing=TRUE), amr.ranges=ramr.tp.nonunique)
  )
  
  RUnit::checkEquals(
    length(plotAMR(data.ranges=ramr.data, amr.ranges=ramr.tp.nonunique)),
    length(GenomicRanges::reduce(ramr.tp.nonunique))
  )
  RUnit::checkEquals(
    length(plotAMR(data.ranges=ramr.data, amr.ranges=ramr.tp.unique)),
    length(GenomicRanges::reduce(ramr.tp.unique))
  )
}
