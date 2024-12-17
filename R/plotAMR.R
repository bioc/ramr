#' Plot aberrantly methylated regions
#'
#' @description
#' `plotAMR` uses `ggplot2` to visualize aberrantly methylated regions (AMRs)
#' at particular genomic locations.
#'
#' @details
#' For every non-overlapping genomic location from `amr.ranges` object,
#' `plotAMR` plots and outputs a line graph of methylation beta values taken
#' from `data.ranges` for all samples from `data.samples`. Samples bearing
#' significantly different methylation profiles ('sample' column of
#' `amr.ranges` object) are highlighted.
#'
#' @param data.ranges A `GRanges` object with genomic locations and
#' corresponding beta values included as metadata.
#' @param data.samples A character vector with sample names (a subset of
#' metadata column names) to be included in the plot. If `NULL` (the default),
#' then all samples (metadata columns) are included.
#' @param amr.ranges An output of `getAMR` - a `GRanges` object that contain
#' aberrantly methylated regions (AMRs).
#' @param highlight An optional list of samples to highlight. If NULL (the
#' default), will contain sample IDs from the `sample` metadata column of
#' `amr.ranges` object.
#' @param title An optional title for the plot. If NULL (the default), plot
#' title is set to a genomic location of particular AMR.
#' @param labs Optional axis labels for the plot. Default: c("genomic position",
#' "beta value").
#' @param window An optional integer constant to expand genomic ranges of the
#' `amr.ranges` object (the default: 300).
#' @return The output is a list of `ggplot` objects.
#' @seealso \code{\link{getAMR}} for identification of AMRs,
#' \code{\link{getUniverse}} for info on enrichment analysis,
#' \code{\link{simulateAMR}} and \code{\link{simulateData}} for the generation
#' of simulated test data sets, and `ramr` vignettes for the description of
#' usage and sample data.
#' @examples
#'   data(ramr)
#'   plotAMR(ramr.data, ramr.samples, ramr.tp.unique[1])
#'   library(gridExtra)
#'   do.call("grid.arrange",
#'           c(plotAMR(ramr.data, ramr.samples, ramr.tp.nonunique), ncol=2))
#' @export
plotAMR <- function (data.ranges,
                     data.samples=NULL,
                     amr.ranges,
                     highlight=NULL,
                     title=NULL,
                     labs=c("genomic position", "beta value"),
                     window=300)
{
  if (!requireNamespace("ggplot2", quietly=TRUE)) stop("ggplot2 is required for plotting. Please install")
  
  if (is.null(data.samples))
    data.samples <- colnames(GenomicRanges::mcols(data.ranges))
  amr.ranges.reduced  <- GenomicRanges::reduce(amr.ranges, min.gapwidth=window, with.revmap=TRUE)
  amr.ranges.relisted <- BiocGenerics::relist(amr.ranges[unlist(amr.ranges.reduced$revmap)], amr.ranges.reduced$revmap)
  plot.list <- list()

  for (i in seq_along(amr.ranges.relisted)) {
    plot.ranges <- unlist(amr.ranges.relisted[i])
    # revmap.rows <- unique(unlist(plot.ranges$revmap))
    data.hits   <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(data.ranges, plot.ranges, maxgap=window, ignore.strand=TRUE)))
    if (length(data.hits)>0) {
      plot.data <- data.table::as.data.table(data.ranges[data.hits, data.samples])
      plot.data$median <- apply(plot.data[, ..data.samples], 1, median, na.rm=TRUE)

      colorify       <- c("median", if (is.null(highlight)) unique(plot.ranges$sample), highlight)
      plot.data.melt <- data.table::melt.data.table(plot.data, id.vars=c("seqnames","start","end","width","strand"),
                                       variable.name="sample", value.name="beta")
      plot.data.melt[, `:=` (
        size=0,
        alpha=0.5,
        color=factor("lightgrey",levels=c("lightgrey", colorify))
      )]
      
      plot.data.melt[sample %in% colorify, `:=` (alpha=0.9, color=sample)]
      for (j in seq_along(plot.ranges)) {
        plot.data.melt[sample==plot.ranges$sample[j] & start %in% GenomicRanges::start( data.ranges[unlist(plot.ranges[j]$revmap)] ),
                       size:=1]
      }
      
      gene.plot <- ggplot2::ggplot(plot.data.melt, ggplot2::aes(x=start, y=beta, group=sample, color=color, alpha=alpha)) +
        ggplot2::geom_line(linewidth=0.5) +
        ggplot2::geom_point(mapping=ggplot2::aes(size=size)) +
        ggplot2::scale_x_continuous(name=labs[1]) +
        ggplot2::scale_y_continuous(name=labs[2], limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
        ggplot2::scale_color_discrete(name="samples", limits=colorify) +
        ggplot2::scale_alpha_continuous(guide="none") +
        ggplot2::scale_size_identity(guide="none") +
        ggplot2::theme_light() +
        ggplot2::theme(legend.text=ggplot2::element_text(size=8),
                       axis.text.x=ggplot2::element_text(size=8, angle=0),
                       axis.text.y=ggplot2::element_text(size=8)) +
        ggplot2::ggtitle(
          if (is.null(title)) as.character(GenomicRanges::reduce(plot.ranges)) else title
        )

      plot.list[length(plot.list)+1] <- list(gene.plot)
    }
  }

  return(plot.list)
}
