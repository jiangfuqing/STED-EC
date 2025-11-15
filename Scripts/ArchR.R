pkgs <- c("Seurat","Signac","ArchR","tidyverse","Matrix","patchwork","grid","Rsamtools",
          "ggalluvial","clusterProfiler","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene",
          "parallel","BSgenome.Mmusculus.UCSC.mm10","TFBSTools","chromVARmotifs","motifmatchr",
          "magrittr","viridis","igraph", "clustree", "reticulate", "ggrepel")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

# accept args from pbs file
sampleId <- "E185-P-12"
Species <- "mouse"

#threads default to 1 in windows
addArchRThreads(threads = 50)

# set work directory
#setwd("~/MISAR-seq")

if (Species=="mouse") {addArchRGenome("mm10")} else {addArchRGenome("hg38")}

#Download E185-P-12.arrow file from google drive
proj <- readRDS("")


proj@sampleColData@listData[["ArrowFiles"]] <- "~/Data/E185-P-12.arrow"
