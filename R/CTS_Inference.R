rm(list=ls())
library(sRACIPE)
#library(LSD)
#library(gplots)
#library(biomaRt)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(STICCC)
#source("HPCFunctions.R")
#source("R/HPCFunctions.R")
set.seed(123)


# global params
topoName <- "CTS"
runSim <- TRUE     
doPCA <- TRUE
forceSTICCC <- FALSE
nSamples <- 10000
pseudocount <- F

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
topo$Type[which(topo$Type %% 2 == 0)] = 2
topo$Type[which(topo$Type %% 2 == 1)] = 1


# simulate topology if needed
if(runSim) {
  racipe <- simTopo(topo, numModels = nSamples)
  saveRDS(racipe, file.path(outputDir, paste0("simData_",topoName,".Rds")))
} else {
  racipe <- readRDS(file.path(outputDir, paste0("simData_",topoName,".Rds")))
}

# normalize data using built in function, or manually
exprMat <- assay(racipe)[,]

if(pseudocount) {
  racipe_norm <- sracipeNormalize(racipe)
  exprMat_norm <- assay(racipe_norm)[,]
} else {
  exprMat_norm <- log2(exprMat)
  means <- rowMeans(exprMat_norm)
  sds <-  apply(exprMat_norm, 1, sd)
  exprMat_norm <- sweep(exprMat_norm, 1, means, FUN = "-")
  exprMat_norm <- sweep(exprMat_norm, 1, sds, FUN = "/")
  
}




# create SCE object
## TODO: make the lines below into a small wrapper method createVIC()
vic <- SingleCellExperiment(assays = SimpleList(counts=exprMat, normcounts=exprMat_norm))
vic@metadata$experimentName <- paste0(topoName)
vic@metadata$topoName <- topoName
vic@metadata$topo <- topo
vic@metadata$params <- list(sample_radius=0.05, plotScalingFactor=1, gridPlotScalingFactor=1, minNeighbors=5, verbose=T)


# add metadata
metadata <- prepMetadata(exprMat_norm, cluster = T, k = 2)
rownames(metadata) <- metadata$SampleID
colData(vic) <- DataFrame(metadata)
colnames(vic) <- colData(vic)$SampleID


# Plot network
plotNetwork(vic)


# run PCA
if(doPCA) {
  pca <- runPCA(assay(vic, "normcounts"))
  saveRDS(pca, file.path(outputDir,"PCA_res.Rds"))
} else {
  pca <- readRDS(file.path(outputDir,"PCA_res.Rds"))
}

# add PCA to SCE object
reducedDim(vic, "PCA") <- pca$x
vic@metadata$pca_data <- pca[1:4]
vic@metadata$pca_summary <- summary(pca)


# compute grid based on PCA
vic <- computeGrid(vic)

# compute pairwise distance between points
vic <- computeDist(vic)


# compute trajectories
vic@metadata$params$useGenes <- T


vic_fname <- file.path(outputDir, paste0("vic_",topoName,".Rds"))
if(!file.exists(vic_fname) | forceSTICCC) {
  vic <- DCComputeTrajectorySCE_2022(vic, v2=T)
  saveRDS(vic, vic_fname)
} else {
  vic <- readRDS(vic_fname)
}


# invert v2 for interpretability
# Multiply in vectors by -1
colData(vic)$dX_in <- -1 * colData(vic)$dX_in
colData(vic)$dY_in <- -1 * colData(vic)$dY_in

# Plot results
minMagnitude <- 0.001
scalingFactor <- 2
arrowheadSize <- 0.5


vic <- computeGrid(vic, grid.length = 15)  

vic <- DCComputeGridVectors(vic, inVectors = F, combine = T, unitVectors = F, how="avg+")

DCPlotGrid_Grey(sce = vic,
                plotLoadings = F,
                loadingFactor = 3.5,
                plotSuffix = paste0("_may24_avg+"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                minimal = F,
                arrowheadSize = arrowheadSize)


### Plot out - in
vic <- DCComputeGridVectors(vic, inVectors = F, combine = T, unitVectors = F, how="avg-")

DCPlotGrid_Grey(sce = vic,
                plotLoadings = F,
                loadingFactor = 3.5,
                plotSuffix = paste0("_may24_avg-"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                minimal=F,
                arrowheadSize = arrowheadSize)


grid_diff_actual <- vic@metadata$grid.df


# Plot v1 & v2 separately
vic <- DCComputeGridVectors(vic, inVectors = T, combine = F, unitVectors = F)

DCPlotGrid_Grey(sce = vic,
                plotLoadings = F,
                loadingFactor = 3.5,
                plotSuffix = paste0("_may24_v2"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                minimal = F,
                arrowheadSize = arrowheadSize)



vic <- DCComputeGridVectors(vic, inVectors = F, combine = F, unitVectors = F)

DCPlotGrid_Grey(sce = vic,
                plotLoadings = T,
                loadingFactor = 5,
                plotSuffix = paste0("_may24_v1_loadings"),
                minMagnitude = minMagnitude,
                scalingFactor = scalingFactor,
                minimal = F,
                arrowheadSize = arrowheadSize)







