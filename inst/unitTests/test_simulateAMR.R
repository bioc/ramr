test_simulateAMR <- function () {
  data(ramr)
  RUnit::checkException(
    simulateAMR("ramr.data", 2)
  )
  RUnit::checkException(
    simulateAMR(ramr.data, nsamples=2, sample.names=c("a"))
  )
  RUnit::checkEquals(
    length( simulateAMR(ramr.data, nsamples=6, exclude.ranges=ramr.tp.unique) ),
    6
  )
  RUnit::checkEquals(
    length( simulateAMR(ramr.data, nsamples=1000, regions.per.sample=1000, samples.per.region=1000, merge.window=1000) ),
    143**2
  )
  
  unique    <- simulateAMR(ramr.data, nsamples=29, regions.per.sample=2, samples.per.region=1)
  nonunique <- simulateAMR(ramr.data, nsamples=20, exclude.ranges=unique, merge.window=1000, regions.per.sample=5, samples.per.region=3)
  noise     <- simulateAMR(ramr.data, nsamples=100, exclude.ranges=c(unique,nonunique), merge.window=1, min.cpgs=1, max.cpgs=1, regions.per.sample=1000, samples.per.region=1)
  RUnit::checkTrue(
    all(unique %outside% nonunique)
  )
  RUnit::checkTrue(
    all(noise %outside% c(unique,nonunique))
  )
  RUnit::checkEquals(
    sapply(list(unique, nonunique, noise), length),
    c(58, 240, 1500)
  )
}