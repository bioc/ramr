# rare long-range methylation aberrations


simulateData <- function (template.ranges,
                          nsamples,
                          amr.ranges=NULL,
                          sample.names=NULL,
                          min.beta=0.001,
                          max.beta=0.999,
                          cores=max(1,parallel::detectCores()-1))
{
  if (!methods::is(template.ranges,"GRanges"))
    stop("'template.ranges' must be a GRanges object")
  if (!is.null(sample.names) & length(sample.names)!=nsamples)
    stop("'sample.names' length must be equal to 'nsamples'")
  
  if (is.null(sample.names))
    sample.names <- paste0("sample", seq_len(nsamples))
  if (!is.null(amr.ranges)) {
    if(!methods::is(amr.ranges,"GRanges"))
      stop("'amr.ranges' must be a GRanges object")
    if ( is.null(amr.ranges$revmap) | !all(unlist(amr.ranges$revmap) %in% seq_along(template.ranges)) )
      stop("Malformed 'amr.ranges' object: 'revmap' field is missing or contains out-of-range values")
    if ( is.null(amr.ranges$sample) | !all(amr.ranges$sample %in% sample.names) )
      stop("Malformed 'amr.ranges' object: 'sample' field is missing or its elements are absent from 'sample.names'")
    if ( !is.numeric(amr.ranges$dbeta) | !all(na.omit(amr.ranges$dbeta)>=0) | !all(na.omit(amr.ranges$dbeta)<=1) )
      stop("Malformed 'amr.ranges' object: 'dbeta' field is missing or is outside the range c(0,1)")
  }
    
  #####################################################################################
  
  getRandomBeta <- function(data.chunk) {
    chunk.filt <- apply(data.chunk, 1, function(x) {
      x.median    <- stats::median(x, na.rm=TRUE)
      x[is.na(x)] <- x.median
      beta.fit    <- suppressWarnings( EnvStats::ebeta(as.numeric(x)) )
      random.beta <- stats::rbeta(nsamples, beta.fit$parameters[1], beta.fit$parameters[2])
      return(random.beta)
    })
    return(t(chunk.filt))
  }
  
  #####################################################################################
  
  template.betas <- as.matrix(mcols(template.ranges, use.names=FALSE))
  chunks <- split(seq_len(nrow(template.betas)), if (cores>1) cut(seq_len(nrow(template.betas)), cores) else 1)
  
  doParallel::registerDoParallel(cores)
  cl <- parallel::makeCluster(cores)
  random.betas <- foreach (chunk=chunks) %dopar% getRandomBeta(template.betas[chunk,])
  random.betas <- do.call(rbind, random.betas)
  colnames(random.betas) <- sample.names
  parallel::stopCluster(cl)
  
  amr.mcols <- if (is.null(amr.ranges)) data.frame() else GenomicRanges::mcols(amr.ranges)
  for (i in seq_len(nrow(amr.mcols))) {
    revmap <- unlist(amr.mcols$revmap[i])
    dbeta  <- sign(0.5 - mean(random.betas[revmap,], na.omit=TRUE)) * amr.mcols$dbeta[i]
    sample <- amr.mcols$sample[i]
    random.betas[revmap, sample] <- random.betas[revmap, sample] + dbeta
  }
  random.betas[random.betas>max.beta] <- max.beta
  random.betas[random.betas<min.beta] <- min.beta
  
  data.ranges <- granges(template.ranges)
  mcols(data.ranges) <- random.betas
  
  return(data.ranges)
}