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
#' @param amr.ranges An output of `getAMR` - a `GRanges` object that contain
#' aberrantly methylated regions (AMRs).
#' @param data.samples A character vector with sample names (a subset of
#' metadata column names) to be included in the plot. If `NULL` (the default),
#' then all samples (metadata columns) are included.
#' @param window An optional integer constant to expand genomic ranges of the
#' `amr.ranges` object (the default: 300).
#' @param ignore.strand Boolean to ignore strand of AMR region. Default: FALSE.
#' @param highlight An optional list of samples to highlight. If NULL (the
#' default), will contain sample IDs from the `sample` metadata column of
#' `amr.ranges` object.
#' @param title An optional title for the plot. If NULL (the default), plot
#' title is set to a genomic location of particular AMR.
#' @param labs Optional axis labels for the plot. Default: c("genomic position",
#' "beta value").
#' @param transform Optional transformation of y-axis. Default: "identity"
#' (no transformation).
#' @param limits Optional limits of y-axis. When default (NULL), limits
#' are c(NA,1) for `transform=="log10"` and c(0,1) otherwise.
#' @param breaks Optional breaks of y-axis. When default (NULL), breaks
#' are `10**(seq(from=-5, to=0, length.out=6))` for `transform=="log10"` and
#' `seq(from=0, to=1, length.out=6)` otherwise.
#' @param verbose Boolean to report progress and timings (default: TRUE).
#' @return The output is a list of `ggplot` objects.
#' @seealso \code{\link{getAMR}} for identification of AMRs,
#' \code{\link{getUniverse}} for info on enrichment analysis,
#' \code{\link{simulateAMR}} and \code{\link{simulateData}} for the generation
#' of simulated test data sets, and `ramr` vignettes for the description of
#' usage and sample data.
#' @examples
#'   data(ramr)
#'   plotAMR(data.ranges=ramr.data, amr.ranges=ramr.tp.unique[1])
#'   library(gridExtra)
#'   do.call("grid.arrange",
#'           c(plotAMR(data.ranges=ramr.data, amr.ranges=ramr.tp.nonunique), ncol=2))
#' @export
plotAMR <- function (data.ranges,
                     amr.ranges,
                     data.samples=NULL,
                     window=300,
                     ignore.strand=FALSE,
                     highlight=NULL,
                     title=NULL,
                     labs=c("genomic position", "beta value"),
                     transform=c("identity", "log1p", "log10"),
                     limits=NULL,
                     breaks=NULL,
                     verbose=TRUE)
{
  if (!requireNamespace("ggplot2", quietly=TRUE) | exists(x="is.test.environment"))
    stop("ggplot2 is required for plotting. Please install")

  transform <- match.arg(transform)
  if (is.null(limits))
    limits <- if (transform=="log10") c(NA,1) else c(0,1)
  if (is.null(breaks))
    breaks <- if (transform=="log10") 10**seq(from=-5, to=0, length.out=6) else seq(from=0, to=1, length.out=6)
  if (is.null(data.samples))
    data.samples <- colnames(GenomicRanges::mcols(data.ranges))

  amr.ranges.reduced  <- GenomicRanges::reduce(amr.ranges, min.gapwidth=window, with.revmap=TRUE, ignore.strand=ignore.strand)
  amr.ranges.relisted <- BiocGenerics::relist(amr.ranges[unlist(amr.ranges.reduced$revmap)], amr.ranges.reduced$revmap)

  if (verbose) message("Plotting ", length(amr.ranges.reduced), " genomic ranges ", appendLF=FALSE)
  tm <- proc.time()

  data.ranges.dt <- data.table::as.data.table(data.ranges[, data.samples], optional=TRUE)
  plot.list <- list()

  for (i in seq_along(amr.ranges.relisted)) {
    plot.ranges <- unlist(amr.ranges.relisted[i])
    if (!is.null(plot.ranges$dbeta))
      plot.ranges <- plot.ranges[order(plot.ranges$dbeta, decreasing=TRUE)]
    # plot.key <- data.table::as.data.table(GenomicRanges::reduce(plot.ranges, ignore.strand=ignore.strand))
    plot.key <- data.table::as.data.table(GenomicRanges::granges(plot.ranges))
    plot.key[, `:=` (start=start-window, end=end+window)]
    data.table::setkeyv(plot.key, c(if (ignore.strand) NULL else "strand", "seqnames", "start", "end"))

    data.hits <- unique(data.table::foverlaps(data.ranges.dt, plot.key, mult="all", nomatch=NULL, which=TRUE)$xid)
    # old.data.hits <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(data.ranges, plot.ranges, maxgap=window, ignore.strand=ignore.strand)))
    # if (!identical(data.hits, old.data.hits)) stop("different hits")
    revmap.hits <- unique(unlist(plot.ranges$revmap, recursive=TRUE, use.names=FALSE))
    if (!all(revmap.hits %in% data.hits))
      stop("AMRs don't map to 'data.ranges'! Was the 'data.ranges' object modified after AMR search?")
    if (length(data.hits)>0) {
      plot.data <- data.ranges.dt[data.hits]
      plot.data$median <- apply(plot.data[, ..data.samples], 1, stats::median, na.rm=TRUE)

      amr.samples <- na.omit(plot.ranges$sample)
      amr.revmaps <- plot.ranges$revmap

      colorify       <- c("median", if (is.null(highlight)) unique(as.character(amr.samples)), highlight)
      plot.data.melt <- data.table::melt.data.table(plot.data, id.vars=c("seqnames","start","end","width","strand"),
                                                    variable.name="sample", value.name="beta")
      plot.data.melt[, `:=` (
        size=0,
        alpha=0.5,
        color=factor("lightgrey",levels=c("lightgrey", colorify))
      )]

      plot.data.melt[sample %in% colorify, `:=` (alpha=0.9, color=sample)]
      for (j in seq_along(plot.ranges)) {
        plot.data.melt[sample %in% amr.samples[j] & start %in% data.ranges.dt[amr.revmaps[[j]], start],
                       size:=1]
      }

      plot.title <- if (is.null(title)) as.character(GenomicRanges::reduce(plot.ranges, min.gapwidth=window, ignore.strand=ignore.strand)) else title
      gene.plot <- ggplot2::ggplot(plot.data.melt, ggplot2::aes(x=start, y=beta, group=sample, color=color, alpha=alpha)) +
        ggplot2::geom_line(linewidth=0.5) +
        ggplot2::geom_point(mapping=ggplot2::aes(size=size)) +
        ggplot2::scale_x_continuous(name=labs[1]) +
        ggplot2::scale_y_continuous(name=labs[2], limits=limits, breaks=breaks, transform=transform) +
        ggplot2::scale_color_discrete(name="samples", limits=colorify) +
        ggplot2::scale_alpha_continuous(guide="none") +
        ggplot2::scale_size_identity(guide="none") +
        ggplot2::guides(color=ggplot2::guide_legend(ncol=1)) +
        ggplot2::theme_light() +
        ggplot2::theme(legend.text=ggplot2::element_text(size=8),
                       axis.text.x=ggplot2::element_text(size=8, angle=0),
                       axis.text.y=ggplot2::element_text(size=8)) +
        ggplot2::ggtitle(plot.title)

      plot.list <- c(plot.list, stats::setNames(c(list(gene.plot)), plot.title))
    }

    if (verbose) cat(sprintf("%3.0f%%\b\b\b\b", 100*i/length(amr.ranges.relisted)))
  }

  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(plot.list)
}
