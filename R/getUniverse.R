getUniverse <- function (data.ranges,
                         merge.window=300,
                         min.cpgs=7)
{
  if (class(data.ranges)!="GRanges")
    stop("'data.ranges' must be be a GRanges object")

  universe.ranges <- GenomicRanges::reduce(data.ranges, min.gapwidth=merge.window, with.revmap=TRUE)
  if (length(universe.ranges)>0) {
    universe.ranges$ncpg <- unlist(lapply(universe.ranges$revmap, length))
    universe.ranges      <- subset(universe.ranges, `ncpg`>=min.cpgs & `width`>2)
  }
  return(universe.ranges)
}
