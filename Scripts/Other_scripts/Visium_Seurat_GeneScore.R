pkgs <- c("Seurat","Signac","ArchR","dior","hdf5r","data.table","dplyr","ggplot2","EnsDb.Mmusculus.v79","EnsDb.Hsapiens.v75","tidyverse","Matrix","patchwork","grid","Rsamtools",
          "ggalluvial","clusterProfiler","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene",
          "parallel","BSgenome.Mmusculus.UCSC.mm10","TFBSTools","chromVARmotifs","motifmatchr",
          "magrittr","viridis","igraph", "clustree", "reticulate", "ggrepel")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

set.seed(123)

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

if (Species=="mouse") {addArchRGenome("mm10")} else {addArchRGenome("hg38")}

# create ArchR object
input_ATAC <- paste0(sampleId,'_filtered_fragments.tsv.bgz')

ArrowFiles <- createArrowFiles(
  inputFiles = input_ATAC,
  sampleNames = sampleId,
  minTSS = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  force = FALSE,
  TileMatParams = list(tileSize = 5000)
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleId,
  copyArrows = FALSE
)

#get imputed gene_score
gene_score <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix") # GeneExpressionMatrix, GeneScoreMatrix
rownames(gene_score) <- rowData(gene_score)$name
#proj <- addImputeWeights(proj, reducedDims = "LSI_RNA", td=1, ka=2, nRep=1, k=3) #Default:td=3,ka=4,nRep=2,k=15.
#gene_score <- imputeMatrix(assay(gene_score), getImputeWeights(proj))

#change gene_scores rownames and colnames
gene_scores <- data.frame(gene_score@assays@data@listData$GeneScoreMatrix)
rownames(gene_scores) <- rowData(gene_score)$name
colnames(gene_scores) <- str_split_fixed(colnames(gene_score),"#",2)[,2]

#add GeneScore assay to dat
dat[["GeneScore"]] <- CreateAssayObject(counts = gene_scores)

#read image file
image <- Read10X_Image(image.dir = "./", image.name = paste0(sampleId,"-HE.png"), filter.matrix = F)

image@coordinates <- image@coordinates[na.omit(match(colnames(dat),  rownames(image@coordinates))),]

rownames(image@coordinates) <- colnames(dat)

image@coordinates$imagerow <- image@coordinates$imagerow - 95
image@coordinates$imagecol <- image@coordinates$imagecol - 95

image <- image[Cells(x = dat)]

DefaultAssay(object = image) <- "RNA"
dat[[sampleId]] <- image

#read RNA counts from h5 file
input_RNA <- 'raw_feature_bc_matrix.h5'

seRNA <- import10xFeatureMatrix(
  input = input_RNA,
  names = sampleId,
  featureType = "Gene Expression" # or "Peaks" for ATAC
)

# add RNA data to proj
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA)

#add "Gex_MitoRatio" and "Gex_RiboRatio" to proj
nUMI <- Matrix::colSums(assay(seRNA))
if(Species=="human"){
MitoRatio <- Matrix::colSums(assay(seRNA)[grep("^Mt", rownames(assay(seRNA))), ])/nUMI
RiboRatio <- Matrix::colSums(assay(seRNA)[grep("^Rp", rownames(assay(seRNA))), ])/nUMI} else{MitoRatio <- Matrix::colSums(assay(seRNA)[grep("^Mt", rownames(assay(seRNA))), ])/nUMI
             RiboRatio <- Matrix::colSums(assay(seRNA)[grep("^Rp", rownames(assay(seRNA))), ])/nUMI}
qcInfo <- DataFrame(MitoRatio = MitoRatio, RiboRatio = RiboRatio)
colnames(qcInfo) <- paste0("Gex_", colnames(qcInfo))

# filter non-tissue barcode
barcode <- barcode[z$barcode,]

# read and filter per_barcode_metrics.csv file which from cellranger
metrics <- read.csv("per_barcode_metrics.csv", header = T, row.names = 1)
metrics <- merge(metrics,barcode,by="row.names")
rownames(metrics) <- metrics$Row.names
rownames(metrics) <- paste0(sampleId,"#",rownames(metrics))
metrics <- merge(metrics,qcInfo,by="row.names")
rownames(metrics) <- metrics$Row.names

metrics <- transform(metrics, percent_mito_atac=(metrics$atac_mitochondrial_reads/metrics$atac_raw_reads)*100,
                     Gex_MitoRatio=Gex_MitoRatio*100, Gex_RiboRatio=Gex_RiboRatio*100, 
                     FRiP=(metrics$atac_peak_region_fragments/metrics$atac_fragments)*100, 
                     percent_TSS_fragments=(metrics$atac_TSS_fragments/metrics$atac_fragments)*100, 
                     total_fragments=metrics$atac_fragments, Sample=sampleId)

QC_index <- c("array_col","array_row","Gex_MitoRatio","Gex_RiboRatio","percent_mito_atac","FRiP","percent_TSS_fragments","total_fragments")

# add QC_index into proj
for (i in QC_index) {
  proj <- addCellColData(ArchRProj = proj, data = metrics[,i], cells = rownames(metrics),name = i, force = T)
}

ColData <- getCellColData(proj)

ColData <- ColData[,c("Sample", "array_col", "array_row","Gex_nUMI","Gex_nGenes","Gex_MitoRatio","Gex_RiboRatio","TSSEnrichment","nFrags","log10_unique_fragments","total_fragments","FRiP","percent_mito_atac","percent_TSS_fragments")]

colnames(ColData) <- str_split_fixed(colnames(ColData),"#",2)[,2]

dat <- AddMetaData(object = dat, metadata = ColData)

DefaultAssay(dat) <- "RNA"

if(Species == "mouse"){
  dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-")
  dat[["percent.rp"]] <- PercentageFeatureSet(dat, pattern = "^Rp[sl]")} else{dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
  dat[["percent.rp"]] <- PercentageFeatureSet(dat, pattern = "^RP[SL]") }

saveRDS(dat, file = paste0(sampleId,"_Visium_Seurat_GeneScore1.rds"))

#save RNA object
names(dat@assays) <- c("Spatial","ATAC","GeneScore")
write_h5(dat,paste0(sampleId,"-RNA.h5"),assay.name="spatial")

#save ATAC object
names(dat@assays) <- c("RNA","Spatial","GeneScore")
write_h5(dat,paste0(sampleId,"-ATAC.h5"),assay.name="spatial")

#save GeneScore object
names(dat@assays) <- c("RNA","ATAC","Spatial")
write_h5(dat,paste0(sampleId,"-GeneScore.h5"),assay.name="spatial")