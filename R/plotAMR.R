#' Plot aberrantly methylated regions
#'
#' @description
#' `plotAMR` uses `ggplot2` to visualize aberrantly methylated regions (AMRs)
#' at particular genomic locations.
#'
#' @details
#' For every non-overapping genomic location from `ramr.ranges` object, `plotAMR`
#' plots and outputs a line graph of methylation beta values taken from `data.ranges`
#' for all samples from `data.samples`. Samples bearing significantly different
#' methylation profiles ('sample' column of `ramr.ranges` object) are highlighted.
#'
#' @param data.ranges A `GRanges` object with genomic locations and corresponding
#' beta values included as metadata.
#' @param data.samples A character vector with sample names (a subset of metadata
#' column names) to be included in the plot.
#' @param ramr.ranges An output of `getAMR` - a `GRanges` object that contain
#' aberrantly methylated regions (AMRs).
#' @param highlight An optional list of samples to highlight. If NULL (default),
#' will contain sample IDs from the `sample` metadata column of `ramr.ranges` object.
#' @param title An optional title for the plot. If NULL (default), plot title is
#' set to a genomic location of particular AMR.
#' @param window An optional integer constant to expand genomic ranges of the
#' `ramr.ranges` object.
#' @return The output is a list of `ggplot` objects.
#' @seealso \code{\link{getAMR}} for identification of AMRs, \code{\link{getUniverse}} for
#' info on enrichment analysis
#' @examples
#' \dontrun{
#'   plotAMR(ramr.data, ramr.samples, ramr.tp.unique[1])
#'
#'   # library(gridExtra)
#'   do.call("grid.arrange", c(plotAMR(ramr.data, ramr.samples, ramr.tp.nonunique), ncol=2))
#' }
#' @import GenomicRanges
#' @import ggplot2
#' @importFrom matrixStats rowMedians
#' @importFrom reshape2 melt

#' @export
plotAMR <- function (data.ranges,
                     data.samples,
                     ramr.ranges,
                     highlight=NULL,
                     title=NULL,
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
      if (is.null(title))
        title <- as.character(reduce(plot.ranges))
      colorify       <- c("median",highlight)
      plot.data.melt <- reshape2::melt(plot.data, id.vars=c("seqnames","start","end","width","strand"),
                                       variable.name="sample", value.name="beta")
      plot.data.melt <- cbind(plot.data.melt, list(alpha=0.5,color=factor("lightgrey",levels=c("lightgrey",colorify))))

      plot.data.melt[plot.data.melt$sample %in% colorify, "alpha"] <- 0.9
      for (sample.name in colorify)
        plot.data.melt[plot.data.melt$sample==sample.name,"color"] <- sample.name

      gene.plot <- ggplot(plot.data.melt, aes_string(x="start", y="beta", group="sample", color="color", alpha="alpha")) +
        geom_line(size=0.5) +
        geom_point(size=1) +
        scale_x_continuous(name="position") +
        scale_y_continuous(name="beta value", limits=c(-0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
        scale_color_discrete(name="samples", limits=colorify) +
        scale_alpha_continuous(guide="none") +
        theme(legend.text=element_text(size=8), #legend.position="none", legend.title=element_blank(),
              axis.text.x=element_text(size=8, angle=0), #, color=.data.cpgs.colors),
              axis.text.y=element_text(size=8)) +
        ggtitle(title)

      # print(gene.plot)
      plot.list[length(plot.list)+1] <- list(gene.plot)
    }
  }

  return(plot.list)
}
