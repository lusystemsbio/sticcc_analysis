##### GLOBALS & SETUP #####
rm(list=ls())
library(sRACIPE)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(STICCC)
library(MASS)
library(pracma)
set.seed(123)


# global params
topoName <- "CTS"
forceSim <- FALSE     
forcePCA <- FALSE
forceSTICCC <- FALSE
saveNetworkPlot <- FALSE
nSamples <- 10000
pseudocount <- T
numClusters <- 4
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# directory setup
topoDir <- file.path(getwd(),topoName)
outputDir = file.path(topoDir,"data")
plotDir <- file.path(getwd(), topoName, paste0(topoName, "_2024"))
if(!dir.exists(topoDir)) {
  dir.create(topoDir)
}
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}
if(!dir.exists(plotDir)) {
  dir.create(plotDir)
}



# load topology file
topo <- loadTopo(topoName)
topo$Type[which(topo$Type %% 2 == 0)] = 2
topo$Type[which(topo$Type %% 2 == 1)] = 1


##### SIMULATE GRN #####
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

##### PCA & STICCC SETUP #####
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

###### IDENTIFY BASINS ######
# Identify basin positions
pca_df <- as.data.frame(reducedDim(stic,"PCA"))
kde_result <- kde2d(pca_df$PC1, pca_df$PC2, n = 100)
kde_density <- kde_result$z
energy_landscape <- -log(kde_density)

# Identify local minima as basins
local_minima <- matrix(FALSE, nrow = nrow(energy_landscape), ncol = ncol(energy_landscape))
for (i in 2:(nrow(energy_landscape) - 1)) {
  for (j in 2:(ncol(energy_landscape) - 1)) {
    if (energy_landscape[i, j] == min(energy_landscape[(i-1):(i+1), (j-1):(j+1)])) {
      local_minima[i, j] <- TRUE
    }
  }
}

# Convert minima to PCA coordinates
minima_coords <- which(local_minima, arr.ind = TRUE)
minima_x <- kde_result$x[minima_coords[, 1]]
minima_y <- kde_result$y[minima_coords[, 2]]

# Plot KDE and basins
kde_df <- data.frame(expand.grid(PC1 = kde_result$x, PC2 = kde_result$y), Density = as.vector(kde_density))

image <- ggplot(kde_df, aes(x = PC1, y = PC2)) +
  geom_raster(aes(fill = Density)) +
  geom_contour(aes(z = Density), color = "black") +
  geom_point(data = data.frame(x = minima_x, y = minima_y), aes(x, y), color = "red", size = 3) +
  scale_fill_viridis_c() +
  theme_minimal() 

density_basins_fname <- file.path(plotDir,paste0(topoName,"_kde2d_basins.pdf"))
pdf(density_basins_fname, height = 10, width = 10)
print(image)
dev.off()


###### EXPRESSION HEATMAP ######
# Plot gene expression distribution (heatmap by cluster)
cluster_order <- order(colData(stic)$Cluster)
ha_df <- data.frame(Cluster=as.character(colData(stic)$Cluster)[cluster_order])

column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]),
                                                          "5"=unname(cbPalette[5]),"6"=unname(cbPalette[6]))))
# Create the heatmap with annotation
image <- Heatmap(exprMat_norm[,cluster_order], 
                 name = "Expression", 
                 top_annotation = column_annotation,
                 row_names_gp=gpar(fontsize=16),
                 cluster_columns = F)

wt_hmap_fname <- file.path(plotDir,paste0(topoName,"_expression_hmap.pdf"))
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()


# compute grid based on PCA
stic <- computeGrid(stic)

# compute pairwise distance between points
stic <- computeDist(stic)


##### RUN STICCC #####
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

##### PLOT RESULTS #####
# Plot results
minMagnitude <- 0.001
scalingFactor <- 2
arrowheadSize <- 0.5


# PCA with cluster annotations
plotVectors(sce = stic,
            colorVar = "Cluster",
            plotLoadings = F,
            plotSuffix = paste0("_jul24_v1_grey"),
            scalingFactor = scalingFactor,
            plotNoVectors = T
)



# Compute gridpoints
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







