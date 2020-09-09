plotAMR <- function (data.ranges,
                      data.samples,
                      ramr.ranges,
                      highlight=NULL,
                      window=300)
{
  ramr.ranges.reduced  <- GenomicRanges::reduce(ramr.ranges, min.gapwidth=window, with.revmap=TRUE)
  ramr.ranges.relisted <- relist(ramr.ranges[unlist(ramr.ranges.reduced$revmap)], ramr.ranges.reduced$revmap)
  plot.list <- list()

  for (i in 1:length(ramr.ranges.relisted)) {
    plot.ranges <- unlist(ramr.ranges.relisted[i])
    revmap.rows <- unique(unlist(plot.ranges$revmap))
    data.hits   <- unique(queryHits(findOverlaps(data.ranges, plot.ranges, maxgap=window, ignore.strand=TRUE)))
    if (length(data.hits)>0) {
      plot.data <- data.frame(data.ranges[data.hits, data.samples], check.names=FALSE, stringsAsFactors=FALSE)
      colnames(plot.data) <- c(colnames(plot.data)[1:5], data.samples)
      plot.data$median <- matrixStats::rowMedians(as.matrix(plot.data[,data.samples]), na.rm=TRUE)

      if (is.null(highlight))
        highlight <- unique(plot.ranges$sample)
      colorify       <- c("median",highlight)
      plot.data.melt <- melt(plot.data, id.vars=c("seqnames","start","end","width","strand"),
                             variable.name="sample", value.name="beta")
      plot.data.melt <- cbind(plot.data.melt, list(alpha=0.5,color=factor("lightgrey",levels=c("lightgrey",colorify))))

      plot.data.melt[plot.data.melt$sample %in% colorify, "alpha"] <- 0.9
      for (sample.name in colorify)
        plot.data.melt[plot.data.melt$sample==sample.name,"color"] <- sample.name

      gene.plot <- ggplot(plot.data.melt, aes(x=start, y=beta, group=sample, color=color, alpha=alpha)) +
        geom_line(size=0.5) +
        geom_point(size=1) +
        scale_x_continuous(name="position") +
        scale_y_continuous(name="beta value", limits=c(-0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
        scale_color_discrete(name="samples", limits=colorify) +
        scale_alpha_continuous(guide="none") +
        theme(legend.text=element_text(size=8), #legend.position="none", legend.title=element_blank(),
              axis.text.x=element_text(size=8, angle=0), #, color=.data.cpgs.colors),
              axis.text.y=element_text(size=8)) +
        ggtitle(as.character(reduce(plot.ranges)))

      # print(gene.plot)
      plot.list[length(plot.list)+1] <- list(gene.plot)
    }
  }

  return(plot.list)
}


# plotRAMR <- function (data.ranges,
#                       data.samples,
#                       ramr.ranges,
#                       window=300,
#                       title=NULL)
# {
#   ramr.ranges.reduced  <- reduce(ramr.ranges, with.revmap=TRUE)
#   ramr.ranges.relisted <- relist(ramr.ranges[unlist(ramr.ranges.reduced$revmap)], ramr.ranges.reduced$revmap)
#
#   for (i in 1:length(ramr.ranges.relisted)) {
#     plot.ranges <- unlist(ramr.ranges.relisted[i])
#     revmap.rows <- unique(unlist(plot.ranges$revmap))
#     data.hits   <- unique(queryHits(findOverlaps(data.ranges, plot.ranges, maxgap=window, ignore.strand=TRUE)))
#     plot.data   <- data.frame(data.ranges[data.hits, data.samples], check.names=FALSE, stringsAsFactors=FALSE)
#     colnames(plot.data) <- c(colnames(plot.data)[1:5], data.samples)
#     plot.data$median <- matrixStats::rowMedians(as.matrix(plot.data[,data.samples]), na.rm=TRUE)
#
#     colorify       <- c("median",unique(plot.ranges$sample))
#     plot.data.melt <- melt(plot.data, id.vars=c("seqnames","start","end","width","strand"),
#                            variable.name="sample", value.name="beta")
#     plot.data.melt <- cbind(plot.data.melt, list(alpha=0.5,color=factor("lightgrey",levels=c("lightgrey",colorify))))
#
#     plot.data.melt[plot.data.melt$sample %in% colorify, "alpha"] <- 0.9
#     for (sample.name in colorify)
#       plot.data.melt[plot.data.melt$sample==sample.name,"color"] <- sample.name
#
#     gene.plot <- ggplot(plot.data.melt, aes(x=start, y=beta, group=sample, color=color, alpha=alpha)) +
#       geom_line(size=0.5) +
#       geom_point(size=1) +
#       # scale_x_discrete(limits=revmap.rows) +
#       scale_x_continuous(name="position") +
#       scale_y_continuous(name="beta value", limits=c(-0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
#       scale_color_discrete(name="samples", limits=colorify) +
#       scale_alpha_continuous(guide="none") +
#       theme(legend.text=element_text(size=8), #legend.position="none", legend.title=element_blank(),
#             axis.text.x=element_text(size=8, angle=0), #, color=.data.cpgs.colors),
#             axis.text.y=element_text(size=8)) +
#       ggtitle( ifelse( !is.null(title), title, as.character(reduce(plot.ranges)) ) )
#
#     print(gene.plot)
#   }
# }
