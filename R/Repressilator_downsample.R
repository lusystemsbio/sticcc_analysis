rm(list = ls())
library(sRACIPE)
library(ggpubr)
library(STICCC)
#source("R/HPCFunctions.R")
#source("R/functions_radius_optimization.R")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# global params
topoName <- "repressilator"

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
fname = file.path(outputDir,paste0("stic_",topoName,".Rds"))
stic_old <- readRDS(fname)
stic_old <- computeGridVectors(stic_old)
grid.df <- stic_old@metadata$grid.df

# plot with full data
plotVectors(stic_old, plotSuffix = "10ksamples_aug2024", scalingFactor = 0.3)


# plot with limited data
plotVectors(stic_old[,sample(1:10000, 1000)], plotSuffix = "1ksamples_aug2024", scalingFactor = 0.3)



# downsample
sample <- sample(1:10000, 200)
stic <- stic_old[,sample]


# compute new trajectories on limited data
stic@metadata$params$sample_radius <- 0.1
stic <- runSTICCC(stic, v2=T)

# plot
stic <- computeGrid(stic)
stic <- computeGridVectors(stic)

stic@metadata$params$xMin <- -3
stic@metadata$params$xMax <- 3
stic@metadata$params$yMin <- -3
stic@metadata$params$yMax <- 3

plotGrid(stic, plotSuffix = "200samples_aug2024")



