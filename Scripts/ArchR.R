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

#Download E185-P-12-SpatialGlue.rds E185-P-12.arrow file from google drive
proj <- readRDS("E185-P-12-SpatialGlue.rds")
proj@sampleColData@listData[["ArrowFiles"]] <- "~/Data/E185-P-12.arrow"

#plot SpatialGlue merged clusters
p1 <- plotEmbedding(proj, name = "ATAC_merge", embedding = "UMAP_ATAC", rastr = FALSE, plotAs = "points", pal = color_list, size = 0.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "RNA_merge", embedding = "UMAP_RNA", rastr = FALSE, plotAs = "points", pal = color_list, size = 0.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "Glue_merge", embedding = "UMAP_Combined", rastr = FALSE, plotAs = "points", pal = color_list, size = 0.5, labelAsFactors=F, labelMeans=F)

plotPDF(p1,p2,p3, name = "UMAP-ATAC-RNA-Combined-SpatialGlue-Merge", addDOC = FALSE)

#get each sample cellcoldata to plot Spatial cluster
idxSample <- list()
cellsSample <- list()
df <- list()
for(i in sampleId){
idxSample[[i]] <- BiocGenerics::which(proj$Sample %in% i)
cellsSample[[i]] <- proj$cellNames[idxSample[[i]]]
df[[i]] <- getCellColData(proj[cellsSample[[i]],])
}

pdf(file = paste0(sampleId,"_ArchR_ATAC-RNA-Combined_Merge_Cluster.pdf"), width=8.6, height=8.6)
for(i in sampleId){
imported_raster=OpenImageR::readImage(paste0(i, "-HE.jpg"))
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
j <- df[[i]][,"ATAC_merge"]
p4 <- ggplot(as.data.frame(df[[i]]), aes(x = as.numeric(array_col), y = as.numeric(array_row), color=j)) +
    scale_color_manual(values = color_list[sort(as.numeric(unique(j)))],
                       labels=paste0("A",sort(as.numeric(unique(j))))) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5, shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"),
      legend.text=element_text(size=20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())
  plot(p4)

j <- df[[i]][,"RNA_merge"]
p5 <- ggplot(as.data.frame(df[[i]]), aes(x = as.numeric(array_col), y = as.numeric(array_row), color=j)) +
    scale_color_manual(values = color_list[sort(as.numeric(unique(j)))],
                       labels=paste0("R",sort(as.numeric(unique(j))))) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5, shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"),
      legend.text=element_text(size=20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())
  plot(p5)
    
j <- df[[i]][,"Glue_merge"]
p6 <- ggplot(as.data.frame(df[[i]]), aes(x = as.numeric(array_col), y = as.numeric(array_row), color=j)) +
    scale_color_manual(values = color_list[sort(as.numeric(unique(j)))],
                       labels=paste0("C",sort(as.numeric(unique(j))))) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5, shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"),
      legend.text=element_text(size=20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())
  plot(p6)    

}
dev.off()

#Plot Spatial marker gene Gene_Score and Expression
markerGenes <- c("Foxf1","Tmem108","Timp4")

#plot marker gene in UMAP
p1 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP_Combined",
  quantCut = c(0.01, 0.95),
  rastr = FALSE,
  size = 0.3,
  plotAs = "points",
  continuousSet = "blueYellow",
  imputeWeights = getImputeWeights(proj)
)

p2 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
  name = markerGenes, 
  embedding = "UMAP_Combined",
  quantCut = c(0.01, 0.95),
  rastr = FALSE,
  size = 0.3,
  plotAs = "points",
  continuousSet = "solarExtra",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = c(p1, p2), 
        name = "Plot-UMAP-Marker-Genes-Combined.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

#Plot Spatial marker gene Gene Score (p1) and Expression (p2)
for (j in markerGenes) {
  pdf(file = paste0(j,"_Spatial_Marker.pdf"), width=8.6, height=8.6)
  for (i in sampleId) {
    barcode <- read.csv("~/Other_scripts/barcode.csv", header = T, row.names = 1)
    rownames(barcode) <- paste0(i, "#", rownames(barcode))
    p.data <- as.data.frame(p1[[j]]$data$color)
    colnames(p.data) <- j
    row.names(p.data) <- as.data.frame(proj@cellColData@rownames)[rownames(p1[[j]]$data),]
    
    data <- merge(p.data,barcode,by="row.names")
    
    imported_raster=OpenImageR::readImage(paste0(i,"-HE.jpg"))
    g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
    p3 <- ggplot(data, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=data[,j])) +
      ggtitle(paste0(i,"_GeneScore_",j)) + scale_color_gradientn(colours = ArchRPalettes$blueYellow) +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      geom_point(size = 2.5, shape = 15) +
      expand_limits(x = 0, y = 0) +
      scale_x_continuous(name = NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
      scale_y_reverse(name = NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
      coord_equal(xlim = c(0,51),ylim = c(51,1)) +
      theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
            axis.text=element_text(size = 0),
            axis.ticks = element_blank(),
            axis.title=element_text(size = 0,face = "bold"),
            legend.text=element_text(size = 20),
            legend.title = element_blank(),
            #legend.title = element_text(colour = "black", size = 15, face="bold"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    plot(p3)
    
    p.data <- as.data.frame(p2[[j]]$data$color)
    colnames(p.data) <- j
    row.names(p.data) <- as.data.frame(proj@cellColData@rownames)[rownames(p2[[j]]$data),]
    data <- merge(p.data,barcode,by="row.names")
    
    p4 <- ggplot(data, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=data[,j])) +
      ggtitle(paste0(i,"_Expression_",j)) + scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      geom_point(size = 2.5, shape = 15) +
      expand_limits(x = 0, y = 0) +
      scale_x_continuous(name = NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
      scale_y_reverse(name = NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
      coord_equal(xlim = c(0,51),ylim = c(51,1)) +
      theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
            axis.text=element_text(size = 0),
            axis.ticks = element_blank(),
            axis.title=element_text(size = 0,face = "bold"),
            legend.text=element_text(size = 20),
            legend.title = element_blank(),
            #legend.title = element_text(colour = "black", size = 15, face="bold"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    plot(p4)
  }
  dev.off()
}
dev.off()
