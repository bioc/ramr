ramr
========

# Introduction

*`ramr`* is an R package for detection of low-frequency aberrant methylation events in large datasets
obtained by methylation profiling using array or high-throughput bisulfite sequencing. In addition, package provides
functions to visualize found aberrantly methylated regions (AMRs), and to generate sets of all possible regions to be used
as reference sets for enrichment analysis.

This readme contains condensed info on *`ramr`* usage. For more, please check function-specific help pages and vignettes.

## Current Features

 * Identification of aberrantly methylated regions (AMRs)
 * AMR visualization
 * Generation of reference sets for third-party analyses (e.g. enrichment)


-------

## Installation

### install via Bioconductor - NOT AVAILABLE YET
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ramr")
```

### Install the latest version via install_github
```r
library(devtools)
install_github("BBCG/ramr", build_vignettes=FALSE,
  repos=BiocManager::repositories(),
  dependencies=TRUE, type="source")
```


-------

### Citing the *`ramr`* package
[Nikolaienko et al., 2020 bioRxiv](https://doi.org/10.1101/2020.12.01.403501)

### The data underlying *`ramr`* manuscript
[Replication Data for: "ramr: an R package for detection of rare aberrantly methylated regions"](https://doi.org/10.18710/ED8HSD)

-------

# How to Use

*`ramr`* methods operate on objects of the class *`GRanges`*. The input object for AMR search must in addition contain metadata columns with sample beta values. A typical input object looks like this:

```
GRanges object with 383788 ranges and 845 metadata columns:
             seqnames    ranges strand |         GSM1235534         GSM1235535         GSM1235536 ...
                <Rle> <IRanges>  <Rle> |          <numeric>          <numeric>          <numeric> ...
  cg13869341     chr1     15865      * |  0.801634776091808  0.846486905008704   0.86732154737116 ...
  cg24669183     chr1    534242      * |  0.834138820071765  0.861974610731835  0.832557979806823 ...
  cg15560884     chr1    710097      * |  0.711275180750356   0.70461945838556  0.699487225634589 ...
  cg01014490     chr1    714177      * | 0.0769098196182058 0.0569443780518647 0.0623154673389864 ...
  cg17505339     chr1    720865      * |  0.876413362222415  0.885593263385521  0.877944732153869 ...
         ...      ...       ...    ... .                ...                ...                ... ...
  cg05615487    chr22  51176407      * |   0.84904178467798  0.836538383875097   0.81568519870099 ...
  cg22122449    chr22  51176711      * |  0.882444486059592  0.870804215405886  0.859269224277308 ...
  cg08423507    chr22  51177982      * |  0.886406345093286  0.882430879852752  0.887241923657461 ...
  cg19565306    chr22  51222011      * | 0.0719084295670266 0.0845209871264646 0.0689074604483659 ...
  cg09226288    chr22  51225561      * |  0.724145303755024  0.696281176451351  0.711459675603635 ...

```

This code shows how to do basic analysis with *`ramr`* using provided data files:

```r
library(ramr)
data(ramr)

# search for AMRs
amrs <- getAMR(ramr.data, ramr.samples, ramr.method="beta", min.cpgs=5,
               merge.window=1000, qval.cutoff=1e-3)

# inspect
amrs
plotAMR(ramr.data, ramr.samples, amrs[1])

# generate the set of all possible genomic regions using sample dataset and
# the same parameters as for AMR search
universe <- getUniverse(ramr.data, min.cpgs=5, merge.window=1000)

# enrichment analysis of AMRs using R library LOLA
library(LOLA)
hg19.coredb <- loadRegionDB(system.file("LOLACore", "hg19", package="LOLA"))
core.hits   <- runLOLA(amrs, universe, hg19.coredb, cores=1, redefineUserSets=TRUE)
```



The input object may be obtained using data from various sources. Here we provide two examples:

### Using data from NCBI GEO

The following code pulls (NB: very large) raw files from NCBI GEO database, performes normalization and creates *`GRanges`* object for further analysis using *`ramr`* (system requirements: 22GB of disk space, 64GB of RAM)
```r
library(minfi)
library(GEOquery)
library(GenomicRanges)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# destination for temporary files
dest.dir <- tempdir()

# downloading and unpacking raw IDAT files
suppl.files <- getGEOSuppFiles("GSE51032", baseDir=dest.dir, makeDirectory=FALSE, filter_regex="RAW")
untar(rownames(suppl.files), exdir=dest.dir, verbose=TRUE)
idat.files  <- list.files(dest.dir, pattern="idat.gz$", full.names=TRUE)
sapply(idat.files, gunzip, overwrite=TRUE)

# reading IDAT files
geo.idat <- read.metharray.exp(dest.dir)
colnames(geo.idat) <- gsub("(GSM\\d+).*", "\\1", colnames(geo.idat))

# processing raw data
genomic.ratio.set <- preprocessQuantile(geo.idat, mergeManifest=TRUE, fixOutliers=TRUE)

# creating the GRanges object with beta values
data.ranges <- granges(genomic.ratio.set)
data.betas  <- getBeta(genomic.ratio.set)
sample.ids  <- colnames(geo.idat)
mcols(data.ranges) <- data.betas

# data.ranges and sample.ids objects are now ready for AMR search using ramr
```

### Using Bismark cytosine report files

```r
library(methylKit)
library(GenomicRanges)

# file.list is a user-defined character vector with full file names of Bismark cytosine report files
file.list

# sample.ids is a user-defined character vector holding sample names
sample.ids

# methylation context string, defines if the reads covering both strands will be merged
context <- "CpG"

# fitting beta distribution (filtering using ramr.method "beta" or "wbeta") requires
# that most of the beta values are not equal to 0 or 1
min.beta <- 0.001
max.beta <- 0.999

# reading and uniting methylation values
meth.data.raw <- methRead(as.list(file.list), as.list(sample.ids), assembly="hg19", header=TRUE,
                          context=context, resolution="base", treatment=rep(0,length(sample.ids)),
                          pipeline="bismarkCytosineReport")
meth.data.utd <- unite(meth.data.raw, destrand=isTRUE(context=="CpG"))

# creating the GRanges object with beta values
data.ranges <- GRanges(meth.data.utd)
data.betas  <- percMethylation(meth.data.utd)/100
data.betas[data.betas<min.beta] <- min.beta
data.betas[data.betas>max.beta] <- max.beta
mcols(data.ranges) <- data.betas

# data.ranges and sample.ids objects are now ready for AMR search using ramr
```


License
---------
Artistic License/GPL
