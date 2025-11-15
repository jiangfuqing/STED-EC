
#count to tsv

library(Seurat)
library(SeuratDisk) 
library(hdf5r)
library(data.table)

# accept args from pbs file
args <- commandArgs(TRUE)

sampleId <- args[1]

plan("multicore", workers = 46)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 5 Gb RAM

set.seed(1234)

dft <- function(dataset) {
  colnames <- names(dataset)
  tempdata <- t(dataset)
  rownames <- tempdata[1,]
  tempdata <- as.data.table(tempdata, row.names = 1)
  names(tempdata) <- rownames
  tempdata <- cbind(colnames, tempdata)
  tempdata
  }

#dat <- readRDS(paste0(sampleId,"_Seurat_filter.rds"))
#mat <- dat@assays$RNA@counts
dat <- Read10X_h5("raw_feature_bc_matrix.h5")
mat <- dat$`Gene Expression`

# read spatial barcode file.
barcode <- read.csv("/public/home/fqjiang_gdl/SPARC/barcode_file/SPARC_barcode.csv",header = T,row.names = 1)

pos <- paste0(barcode$array_col, "x", barcode$array_row)
colnames(mat) <- pos
mat <- mat[,match(gtools::mixedsort(colnames(mat)), colnames(mat))]
mat <- tibble::rownames_to_column(as.data.frame(as.matrix(mat)), "GENE")

mat <- dft(mat)
mat <- tibble::column_to_rownames(mat, var = "colnames")

write.table(mat, file = paste0(sampleId,".tsv"), sep="\t" ,quote=F, col.names = F)
