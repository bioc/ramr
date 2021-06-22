test_simulateData <- function () {
  data(ramr)
  RUnit::checkException(
    simulateData("ramr.data", 2)
  )
  RUnit::checkException(
    simulateData(ramr.data, nsamples=2, sample.names=c("a"))
  )
  RUnit::checkException(
    simulateData(ramr.data, nsamples=2, amr.ranges=c("a"))
  )
  RUnit::checkException(
    simulateData(ramr.data[1:100], nsamples=2, amr.ranges=ramr.tp.unique)
  )
  RUnit::checkException(
    simulateData(ramr.data, nsamples=2, amr.ranges=ramr.tp.unique)
  )
  RUnit::checkException(
    simulateData(ramr.data, nsamples=100, amr.ranges=ramr.tp.unique)
  )
  
  noise <- simulateAMR(ramr.data, nsamples=10, merge.window=1, min.cpgs=1, max.cpgs=1,
                       regions.per.sample=100, samples.per.region=1, dbeta=1)
  betas <- as.matrix(GenomicRanges::mcols( simulateData(ramr.data, nsamples=10, amr.ranges=noise, cores=2) ))
  RUnit::checkEquals(
    dim(betas),
    c(3000, 10)
  )
  RUnit::checkTrue(
    sum( betas==0.999 | betas==0.001 ) >= length(noise)
  )
  noise <- simulateAMR(ramr.data, nsamples=10, merge.window=1, min.cpgs=1, max.cpgs=1,
                       regions.per.sample=100, samples.per.region=1, dbeta=NA)
  betas <- as.matrix(GenomicRanges::mcols( simulateData(ramr.data, nsamples=10, amr.ranges=noise, cores=2) ))
  RUnit::checkTrue(
    sum(is.na(betas)) == length(noise)
  )
}