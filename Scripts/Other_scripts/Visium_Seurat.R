
#load packages

library(Seurat)
library(dior) 
library(hdf5r)
library(data.table)
library(dplyr)
library(ggplot2)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v75)
library(Rsamtools)

# accept args from pbs file
args <- commandArgs(TRUE)

sampleId <- args[1]
Species <- args[2]

plan("multicore", workers = 46)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 5 Gb RAM

set.seed(1234)

# read counts
counts <- Read10X_h5("raw_feature_bc_matrix.h5")

# create a Seurat object
dat <- CreateSeuratObject(
  project = sampleId,
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

dat[["ATAC"]] <- CreateAssayObject(counts = counts$Peaks)

#Create assay object from fragments
#fragpath <- CreateFragmentObject("fragments.tsv.gz")

# Get gene annotations for hg38 or mm10
#if(Species == "mouse"){annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)} else{annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)}

# If occur Error:in download.file, please update your GenomeInfoDb version > 1.29.10, or use below scripts
#annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
#atac <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"),fragments = fragpath,annotation = annotation)

#frag_counts = Signac::FeatureMatrix(fragments = fragpath, features = atac@annotation) #ranges

#frag.assay <- Signac::CreateChromatinAssay(counts = frag_counts, min.features = 0, fragments = fragpath)

#dat[["Frags"]] <- CreateAssayObject(counts = frag.assay)

#read image file
image <- Read10X_Image(image.dir = "./", image.name = paste0(sampleId,"-HE.png"), filter.matrix = F)

image@coordinates <- image@coordinates[na.omit(match(colnames(dat),  rownames(image@coordinates))),]

rownames(image@coordinates) <- colnames(dat)

image@coordinates$imagerow <- image@coordinates$imagerow - 95
image@coordinates$imagecol <- image@coordinates$imagecol - 95

image <- image[Cells(x = dat)]

DefaultAssay(object = image) <- "RNA"
dat[[sampleId]] <- image

# read position file for filter
location <- read.table(paste0("position_", sampleId,".txt"), sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])[-1]
x <- as.data.frame(x)
colnames(x) <- "array"
y <- read.csv("/public/home/fqjiang_gdl/SPARC/barcode_file/SPARC_barcode_filter.csv", header = T) 
z <- merge(x,y,by="array")

# read the spatial barcode
barcode <- read.csv("/public/home/fqjiang_gdl/SPARC/barcode_file/SPARC_barcode.csv",header = T,row.names = 1)
#rownames(barcode) <- str_split_fixed(rownames(barcode), "#", 2)[,2]
dat <- AddMetaData(object = dat, metadata = barcode)

# filter non-tissue barcodes
dat <- subset(dat, cells=z$barcode)

# read and filter per_barcode_metrics.csv file which from cellranger
metrics <- read.csv("per_barcode_metrics.csv",header = T,row.names = 1)

metrics <- merge(metrics,barcode,by="row.names")

metrics <- metrics[metrics$Row.names %in% z$barcode,]

rownames(metrics) <- metrics$Row.names

metrics <- transform(metrics, atac_percent_mito=(metrics$atac_mitochondrial_reads/metrics$atac_raw_reads)*100, FRiP=(metrics$atac_peak_region_fragments/metrics$atac_fragments)*100, percent_TSS_fragments=(metrics$atac_TSS_fragments/metrics$atac_fragments)*100, log10_unique_fragments=log10(metrics$atac_fragments), Sample=sampleId)

metrics <- metrics[,c("FRiP", "log10_unique_fragments", "atac_percent_mito", "percent_TSS_fragments", "Sample", "array_col", "array_row" )]

dat <- AddMetaData(object = dat, metadata = metrics)


DefaultAssay(dat) <- "RNA"

# calculate percent of mitochondrial for RNA
if(Species == "mouse"){
  dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-")
  dat[["percent.rp"]] <- PercentageFeatureSet(dat, pattern = "^Rp[sl]")} else{dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
  dat[["percent.rp"]] <- PercentageFeatureSet(dat, pattern = "^RP[SL]") }

#dat <- SCTransform(dat, do.scale = T, verbose = FALSE, assay = "RNA") %>% RunPCA(verbose = FALSE)%>% FindNeighbors(dims = 1:15) %>% FindClusters(resolution = c(0.5,0.8,1.0)) %>% RunUMAP(dims = 1:15, reduction.name = "umap.rna")

pdf(paste0(sampleId,"_Spatial_Feature.pdf"))
options(repr.plot.width = 8, repr.plot.height = 8)

for (i in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "FRiP", "log10_unique_fragments", "atac_percent_mito", "percent_TSS_fragments")){
p <- SpatialPlot(dat, features = i, pt.size.factor = 2, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 1, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)

plot(p)

}

dev.off()

#options(repr.plot.width = 8, repr.plot.height = 8)
#p <- SpatialPlot(dat, label = FALSE, label.size = 3, crop = F, group.by = 'Clusters_RNA_0.8', pt.size.factor = 1.5, image.alpha = 1, stroke = 0)
#p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
#p

saveRDS(dat, file = paste0(sampleId,"_Visium_Seurat.rds"))

#save RNA object
names(dat@assays) <- c("Spatial","ATAC")
write_h5(dat,paste0(sampleId,"-RNA.h5"),assay.name="spatial")

#save ATAC object
names(dat@assays) <- c("RNA","Spatial")
write_h5(dat,paste0(sampleId,"-ATAC.h5"),assay.name="spatial")
