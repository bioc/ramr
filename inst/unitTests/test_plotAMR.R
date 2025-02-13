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

  ### if suggested library is not available
  test.env <- new.env()
  assign(x="is.test.environment", value=TRUE, envir=test.env)
  test.func <- function(f, env, ...) {
    environment(f) <- env
    f(...)
  }
  RUnit::checkException(
    test.func(f=plotAMR, env=test.env, data.ranges=ramr.data, amr.ranges=ramr.tp.nonunique)
  )
}
