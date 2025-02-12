test_simulateData <- function () {
  data(ramr)
  RUnit::checkException(
    simulateData("ramr.data", 2)
  )
  RUnit::checkException(
    simulateData(template.ranges=ramr.data, nsamples=2, sample.names=c("a"))
  )
  RUnit::checkException(
    simulateData(template.ranges=ramr.data, nsamples=2, amr.ranges=c("a"))
  )
  RUnit::checkException(
    simulateData(template.ranges=ramr.data[1:100], nsamples=2, amr.ranges=ramr.tp.unique)
  )
  RUnit::checkException(
    simulateData(template.ranges=ramr.data, nsamples=2, amr.ranges=ramr.tp.unique)
  )
  RUnit::checkException(
    simulateData(template.ranges=ramr.data, nsamples=100, amr.ranges=ramr.tp.unique)
  )
  RUnit::checkException(
    simulateData(template.ranges=ramr.data, nsamples=99, amr.ranges=ramr.tp.unique)
  )
  
  noise <- simulateAMR(template.ranges=ramr.data, nsamples=10, merge.window=1, min.cpgs=1, max.cpgs=1,
                       regions.per.sample=100, samples.per.region=1, dbeta=1)
  betas <- as.matrix(GenomicRanges::mcols( simulateData(template.ranges=ramr.data, nsamples=10, amr.ranges=noise) ))
  RUnit::checkEquals(
    dim(betas),
    c(3000, 10)
  )
  RUnit::checkTrue(
    sum( betas==1 | betas==0 ) >= length(noise)
  )
  noise <- simulateAMR(template.ranges=ramr.data, nsamples=10, merge.window=1, min.cpgs=1, max.cpgs=1,
                       regions.per.sample=100, samples.per.region=1, dbeta=NA)
  betas <- as.matrix(GenomicRanges::mcols( simulateData(ramr.data, nsamples=10, amr.ranges=noise) ))
  RUnit::checkTrue(
    sum(is.na(betas)) == length(noise)
  )
  
  ### tests to cover 0/1/NA values during simulation
  
  data.mcols <- GenomicRanges::mcols(ramr.data)
  set.seed(1)
  data.mcols[sample(x=seq_len(nrow(data.mcols)), size=10), sample(x=seq_len(ncol(data.mcols)), size=10)] <- NA
  data.mcols[sample(x=seq_len(nrow(data.mcols)), size=10), sample(x=seq_len(ncol(data.mcols)), size=10)] <- 0
  data.mcols[sample(x=seq_len(nrow(data.mcols)), size=10), sample(x=seq_len(ncol(data.mcols)), size=10)] <- 1
  data.mcols[1500, ] <- NA
  mod.data <- ramr.data
  GenomicRanges::mcols(mod.data) <- data.mcols
  
  simulateData(mod.data, nsamples=100, compute="beta+binom", compute.estimate="amle", compute.weights="logInvDist")
  RUnit::checkException(
    simulateData(mod.data, nsamples=100, compute="beta+binom", compute.estimate="nmle", compute.weights="logInvDist")
  )
  
}