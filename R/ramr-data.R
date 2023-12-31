#' Simulated Illumina HumanMethylation 450k data set with 3000 CpGs and 100
#' samples
#'
#' Data was simulated using GSE51032 data set as described in the reference.
#' Current data set (\code{"ramr.data"}) contains beta values for 10000 CpGs
#' and 100 samples (\code{"ramr.samples"}), and carries 6 unique
#' (\code{"ramr.tp.unique"}) and 15 non-unique (\code{"ramr.tp.nonunique"})
#' true positive AMRs containing at least 10 CpGs with their beta values
#' increased/decreased by 0.5.
#'
#' @docType data
#'
#' @aliases ramr.data ramr.samples ramr.tp.unique ramr.tp.nonunique
#'
#' @usage data(ramr)
#'
#' @format Objects of class \code{"GRanges"} (\code{"ramr.data, ramr.tp.unique,
#' ramr.tp.nonunique"}) and \code{"character"} (\code{"ramr.samples"}).
#'
#' @keywords data sets
#'
#' @references Nikolaienko et al., 2020
#' (\href{https://doi.org/10.1101/2020.12.01.403501}{bioRxiv})
#'
#' @examples
#'   data(ramr)
#'   amrs <- getAMR(ramr.data, ramr.samples, ramr.method="beta", min.cpgs=5,
#'                  merge.window=1000, qval.cutoff=1e-3, cores=2)
#'   plotAMR(ramr.data, ramr.samples, amrs[1])
#'   plotAMR(ramr.data, ramr.samples, ramr.tp.nonunique[4],
#'           highlight=c("sample7","sample8","sample9"))
#'
"ramr.data" #c("ramr.data","ramr.samples","ramr.tp.unique","ramr.tp.nonunique")
