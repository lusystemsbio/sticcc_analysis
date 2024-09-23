rm(list=ls())
library(sRACIPE)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(STICCC)
set.seed(123)


# global params
topoName <- "repressilator"
forceSim <- FALSE     
forcePCA <- FALSE
forceSTICCC <- FALSE
saveNetworkPlot <- FALSE
nSamples <- 10000
pseudocount <- T
numClusters <- 6

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
sim_fname <- file.path(outputDir, paste0("simData_",topoName,".Rds"))
if(forceSim | !file.exists(sim_fname)) {
  racipe <- simTopo(topo, numModels = nSamples)
  saveRDS(racipe, sim_fname)
} else {
  racipe <- readRDS(sim_fname)
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
stic <- sticSE(topo = topo, exprMat = exprMat, normData = exprMat_norm,
              topoName = topoName, expName = paste0(topoName, "_2024"))


# add metadata
stic <- prepMetadata(stic, exprMat_norm, cluster = T, k = numClusters)


# Plot network
if(saveNetworkPlot) {
  plotNetwork(stic)  
}


# run PCA
stic <- runPCA(stic, save=T, overwrite=forcePCA)

### sanity check
# pcadf <- reducedDim(stic,"PCA")
# ggplot(data=pcadf, aes(x=PC1,y=PC2)) +
#   geom_point()

# compute grid based on PCA
stic <- computeGrid(stic)

# compute pairwise distance between points
stic <- computeDist(stic)



# compute trajectories
stic_fname <- file.path(outputDir, paste0("stic_",topoName,".Rds"))
if(!file.exists(stic_fname) | forceSTICCC) {
  stic <- runSTICCC(stic, v2=T, invertV2=T)
  saveRDS(stic, stic_fname)
} else {
  stic <- readRDS(stic_fname)
}

# invert v2 for interpretability
# Multiply in vectors by -1
#colData(stic)$dX_in <- -1 * colData(stic)$dX_in
#colData(stic)$dY_in <- -1 * colData(stic)$dY_in

# Plot results
minMagnitude <- 0.001
scalingFactor <- 1
arrowheadSize <- 0.5


stic <- computeGrid(stic, grid.length = 15)  



### Plot v1
stic <- computeGridVectors(stic, inVectors = F, combine = F, unitVectors = F, how=NA)

plotGrid(sce = stic,
         colorVar = NA,
         plotLoadings = F,
         plotSuffix = paste0("_jul24_v1_grey"),
         minMagnitude = minMagnitude,
         scalingFactor = scalingFactor,
         arrowheadSize = arrowheadSize
)


# Plot v2
stic <- computeGridVectors(stic, inVectors = T, combine = F, unitVectors = F, how=NA)

plotGrid(sce = stic,
         colorVar = NA,
         plotLoadings = F,
         plotSuffix = paste0("_jul24_v2_grey"),
         minMagnitude = minMagnitude,
         scalingFactor = scalingFactor,
         arrowheadSize = arrowheadSize
)



### Plot net
stic <- computeGridVectors(stic, inVectors = F, combine = T, unitVectors = F, how="net")

plotGrid(sce = stic,
         colorVar = NA,
         plotLoadings = F,
         plotSuffix = paste0("_jul24_net_grey"),
         minMagnitude = minMagnitude,
         scalingFactor = scalingFactor,
         arrowheadSize = arrowheadSize
)



### Plot rev
stic <- computeGridVectors(stic, inVectors = F, combine = T, unitVectors = F, how="rev")

plotGrid(sce = stic,
         colorVar = NA,
         plotLoadings = F,
         plotSuffix = paste0("_jul24_rev_grey"),
         minMagnitude = minMagnitude,
         scalingFactor = scalingFactor,
         arrowheadSize = arrowheadSize
)





### Plot v1 again with loadings
stic <- computeGridVectors(stic, inVectors = F, combine = F, unitVectors = F)

plotGrid(sce = stic,
         colorVar = "Cluster",
         plotLoadings = T,
         plotSuffix = paste0("_jul24_v1_loadings_grey"),
         minMagnitude = minMagnitude,
         scalingFactor = scalingFactor,
         arrowheadSize = arrowheadSize
)










