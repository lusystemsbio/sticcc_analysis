rm(list=ls())
library(sRACIPE)
#library(LSD)
library(gplots)
#library(biomaRt)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(Seurat)

#source("R/HPCFunctions.R")
set.seed(123)
library(STICCC)



# global params
topoName <- "emt_net"
runSim <- FALSE     # whether to simulate topology
doPCA <- FALSE
treatment <- "TGFB1"
timepointList <- c("0d","8h","1d","3d","7d","8h_rm","1d_rm","3d_rm")

# directory setup
topoDir <- file.path(getwd(),topoName)
outputDir = file.path(topoDir,"data")
if(!dir.exists(topoDir)) {
  dir.create(topoDir)
}
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}


# load topology file
topo <- loadTopo(topoName)



# pull in expression data
expr <- read.csv(file.path(outputDir, paste0("GSE147405_A549_",treatment,"_TimeCourse_UMI_matrix.csv")), row.names = 1)
metadata <- read.csv(file.path(outputDir, paste0("GSE147405_A549_",treatment,"_TimeCourse_metadata.csv")), row.names = 1)

# remove non-captured genes from topology
genelist <- unique(c(topo$Source, topo$Target))
genelist_expr <- genelist[which(genelist %in% rownames(expr))]

# rank expressed genes by variance
var_df <- data.frame(row.names = genelist_expr)
vars <- rowVars(as.matrix(expr[genelist_expr,]), useNames = T)
means <- rowMeans(as.matrix(expr[genelist_expr,]))


# manually remove low expression, low variance genes
rm_genes <- c("E2F2","IRF9","PPARD","SRF","KLF4")
genelist_expr <- genelist_expr[which(!genelist_expr %in% rm_genes)]
topo <- topo[which(topo$Source %in% genelist_expr & topo$Target %in% genelist_expr),]



# Use Seurat for PCA on top variable genes
seurat <- CreateSeuratObject(counts = expr)
seurat <- NormalizeData(seurat)

expr_norm <- as.matrix(SeuratObject::LayerData(seurat,"data"))

seurat <- FindVariableFeatures(object = seurat, mean.function = ExpMean, dispersion.function = LogVMR)
seurat <- ScaleData(object = seurat)
seurat <- RunPCA(object = seurat, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
Idents(seurat) <- metadata$Time
PCAPlot(object = seurat)


# create SCE object
## TODO: make the lines below into a small wrapper method createVIC()
stic <- SingleCellExperiment(assays = SimpleList(counts=expr, normcounts=expr_norm))
stic@metadata$experimentName <- paste0(topoName,"_A549_",treatment,"_fwd_globalPCA")
stic@metadata$topoName <- topoName
stic@metadata$topo <- topo
stic@metadata$params <- list(sample_radius=0.2, plotScalingFactor=1, gridPlotScalingFactor=1, minNeighbors=5, verbose=T)


# add metadata
#metadata <- readRDS(file.path(getwd(),topoName,"data","metadata_det.Rds"))
#metadata <- prepMetadata(exprMat_norm, cluster = T, k = 6)
#table(metadata$Cluster)
metadata$SampleID <- rownames(metadata)
colData(stic) <- DataFrame(metadata)
colnames(stic) <- colData(stic)$SampleID


# add PCA to SCE object
reducedDim(vic, "PCA") <- as.data.frame(seurat@reductions$pca@cell.embeddings)
colnames(reducedDim(vic,"PCA")) <- gsub("_", "", colnames(reducedDim(vic,"PCA")))
# 1-4: sdev, rotation, scale, center
#vic@metadata$pca_data <- pca[1:4]
vic@metadata$pca_data <- list(seurat@reductions$pca@stdev, seurat@reductions$pca@feature.loadings)

# summary: needs to have $importance
#vic@metadata$pca_summary <- summary(pca)
impdf <- data.frame(row.names=c("Standard deviation", "Proportion of variance", "Cumulative proportion"))
cumsum <- 0
stdev_sum <- sum(seurat@reductions$pca@stdev)
for (pc in seq_along(colnames(seurat@reductions$pca))) {
  stdev <- seurat@reductions$pca@stdev[pc]
  cumsum <- cumsum + stdev
  impdf[1,paste0("PC",pc)] <- stdev
  impdf[2,paste0("PC",pc)] <- seurat@reductions$pca@stdev[pc] / stdev_sum
  impdf[3,paste0("PC",pc)] <- cumsum / stdev_sum
}
vic@metadata$pca_summary <- list(importance = impdf)



# compute grid based on PCA
vic <- computeGrid(vic)


# remove timepoints after signal removal
vic <- vic[,which(metadata$Time %in% timepointList[1:5])]


# compute pairwise distance between points
vic <- computeDist(vic)


# save SCE object
saveRDS(vic, file.path(outputDir,paste0("vic_",topoName,"_A549_",treatment,"_fwd_globalPCA.Rds")))
















vic <- readRDS(file.path(outputDir,paste0("vic_",topoName,"_A549_",treatment,"_fwd_globalPCA.Rds")))






vic@metadata$experimentName <-  paste0("emt_A549_",treatment,"_fwd_globalPCA_nov22")
vic@metadata$params$nPCs <- 15


# compute grid based on PCA
vic <- computeGrid(vic)
# compute pairwise distance between points
vic@metadata$params$sample_radius <- 0.3
vic@metadata$params$minNeighbors <- 30
vic <- computeDist(vic, numPCs = 15)


# compute out vectors
#undebug(DCComputeTrajectorySCE_2022)
vic <- DCComputeTrajectorySCE_2022(vic, v2=T)


#vic_ss <- DCComputeTrajectorySCE_2022(vic[,sample(3133, 100)], v2=T)

# compute in vectors
#vic <- DCComputeTrajectorySCE_in(vic)

# save to RDS
saveRDS(vic, file = file.path(getwd(),topoName,"data",paste0("vic_",topoName,"_inVectors_fwd_globalPCA_nov22.Rds")))


#### Plot results
minMagnitude <- 0.0001
scalingFactor <- 1500
arrowheadSize <- 0.4


vic <- readRDS(file = file.path(getwd(),topoName,"data",paste0("vic_",topoName,"_inVectors_fwd_globalPCA_nov22.Rds")))
#vic@metadata$params$minNeighbors <- 20

# Multiply in vectors by -1
colData(vic)$dX_in <- -1 * colData(vic)$dX_in
colData(vic)$dY_in <- -1 * colData(vic)$dY_in

colData(vic)$Time <- factor(colData(vic)$Time, levels = timepointList)
colData(vic)$CDH1 <- expr_norm["CDH1", colnames(vic)]
colData(vic)$VIM <- expr_norm["VIM", colnames(vic)]
colData(vic)$JUN <- expr_norm["JUN", colnames(vic)]
colData(vic)$TGFBI <- expr_norm["TGFBI", colnames(vic)]
colorGene <- "TGFBI"



ggplot(pca, aes(x=PC_1, y=PC_2)) +
  geom_point(aes(color=expr_norm["TGFBI",])) +
  scale_color_gradient(low="gray", high="red", name="TGFB1")


# save plots
vic <- computeGrid(vic, 15)
vic <- DCComputeGridVectors(vic, inVectors = F, combine = T, unitVectors = F, how="avg+")
# weird point in a549 TGFB1: Mix2_CGGGTCAGTTCGTGAT
# normal point for comparison: Mix2_AACTCAGTCGAACTGT
#debug(DCPLotLocal)
#DCPLotLocal(vic, "Mix2_AACTCAGTCGAACTGT")

#undebug(DCPlotGrid_Grey)
DCPlotGrid_Grey(sce = vic,
                plotLoadings = F,
                loadingFactor = 3.5,
                plotSuffix = paste0("aug2024_avg+"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                arrowheadSize = arrowheadSize,
                minimal = F)


DCPlotGrid(sce = vic,
           plotLoadings = F,
           loadingFactor = 3.5,
           plotSuffix = paste0("aug2024_avg+"),
           minMagnitude = minMagnitude,
           scalingFactor = scalingFactor,
           arrowheadSize = arrowheadSize,
           minimal = F,
           colorVar = colorGene,
           pointSize = 2)

### Plot out - in
vic <- DCComputeGridVectors(vic, inVectors = F, combine = T, unitVectors = F, how="avg-")

DCPlotGrid_Grey(sce = vic,
                plotLoadings = F,
                loadingFactor = 3.5,
                plotSuffix = paste0("aug2024_avg-"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                arrowheadSize = arrowheadSize,
                minimal = F)



DCPlotGrid(sce = vic,
           plotLoadings = F,
           loadingFactor = 3.5,
           plotSuffix = paste0("aug2024_avg-"),
           minMagnitude = minMagnitude,
           scalingFactor = scalingFactor,
           arrowheadSize = arrowheadSize,
           minimal = F,
           colorVar = colorGene,
           pointSize = 2)


# Plot v1 & v2 separately
vic <- DCComputeGridVectors(vic, inVectors = T, combine = F, unitVectors = F)

DCPlotGrid_Grey(sce = vic,
                plotLoadings = F,
                loadingFactor = 3.5,
                plotSuffix = paste0("nov22_v2"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                arrowheadSize = arrowheadSize,
                minimal = F)


vic <- DCComputeGridVectors(vic, inVectors = F, combine = F, unitVectors = F)

DCPlotGrid_Grey(sce = vic,
                plotLoadings = F,
                loadingFactor = 3.5,
                plotSuffix = paste0("nov22_v1"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                arrowheadSize = arrowheadSize,
                minimal = F)





# Also plot colored by timepoint
vic <- DCComputeGridVectors(vic, inVectors = F, combine = T, unitVectors = F, how="avg-")

DCPlotGrid(sce = vic,
           plotLoadings = F,
           loadingFactor = 3.5,
           plotSuffix = paste0("aug22_avg-_byTime"),
           minMagnitude = minMagnitude,
           scalingFactor = scalingFactor,
           colorVar = "Time")




