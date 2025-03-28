####### Setup ####### 

rm(list=ls())
library(sRACIPE)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(STICCC)
set.seed(123)

# global params
topoName <- "CTS"
set.seed(123)
pc1 <- 1
pc2 <- 2
simTime <- 80
printStart <- 0
printInterval <- 0.05
useIC <- TRUE
clust <- c(2)
iNoise <- 0.1
fc <- 50
doTimeSeries <- F
forceSTICCC <- F
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

# read in simulation results
racipe <- readRDS(file.path(outputDir, paste0("simData_",topoName,".Rds")))


# normalize data using built in function
racipe_norm <- sracipeNormalize(racipe)
genelist <- unique(rownames(racipe_norm))
exprMat <- assay(racipe)
exprMat_norm <- assay(racipe_norm)
simExp <- assay(racipe, 1)
simExp <- log2(1+simExp)
tmpMeans <- rowMeans(simExp)
tmpSds <- apply(simExp,1,sd)

# Read stic object & PCA in
stic <- readRDS(file.path(outputDir, paste0("stic_",topoName,".Rds")))
pca <- readRDS(file.path(outputDir, "PCA_res.Rds"))



####### Simulate transitions ####### 

## Select models from one cluster
selected_clusters <- clust
clust4models <- which(colData(stic)$Cluster %in% selected_clusters)
simData_clust4 <- as.data.frame(t(assay(racipe_norm[,clust4models])))
#summary(simData_clust4$A)
#summary(simData_clust4$B)
#summary(simData_clust4$C)
#summary(simData_clust4$D)
## Cluster 4 has high expression of A, so we will increase production of D


## Create new RACIPE object with initial conditions as cluster 4
subRacipe <- RacipeSE(racipe)
numModels <- dim(subRacipe)[2]



## zero-time simulation to gen params and ics
subRacipe <- sracipeSimulate(subRacipe, numModels = numModels,
                             plots = FALSE, genIC = FALSE,
                             genParams = TRUE, integrate = TRUE,
                             integrateStepSize = 0.01,
                             simulationTime = 0, printInterval = 0.05, printStart = 0, 
                             initialNoise = iNoise, nNoise = 1, simDet = F, anneal = F
)



## Set initial conditions and parameters
sracipeIC(subRacipe) <- assay(racipe,1)
sracipeParams(subRacipe) <- sracipeParams(racipe)

## Change production rate of gene D 
sracipeConfig(subRacipe)$simParams["sigGene"] <- 3
sracipeConfig(subRacipe)$simParams["maxFC"] <- fc




### helper functions
processDataTS <- function(racipe,model,pca,tmpMeans,tmpSds) {
  tsData_model <- t(racipe@metadata$timeSeries)
  tsData_model <- data.frame(tsData_model, Time=row.names(tsData_model), Model=model)
  nCols <- ncol(tsData_model)-2
  
  tsData_model[,1:nCols] <- apply(tsData_model[,1:nCols], 2, listToNumeric)
  
  tsData_model[,1:nCols] <- log2(tsData_model[,1:nCols]+1) # Log transform
  tsData_model[,1:nCols] <- sweep(tsData_model[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  tsData_model[,1:nCols] <- sweep(tsData_model[,1:nCols], 2, tmpSds, FUN = "/") # scale
  newPca <- (scale(tsData_model[,1:nCols], pca$center, pca$scale) %*%
               pca$rotation)
  
  tsData_model$PC1 <- newPca[,1]
  tsData_model$PC2 <- newPca[,2]
  tsData_model$PC3 <- newPca[,3]
  
  return(tsData_model)
}

listToNumeric <- function(x) {
  return(as.numeric(unname(unlist(x))))
}



## Time series simulations of each model separately
tsData <- data.frame()
#sink("/dev/null")
models <- clust4models
i <- 1
if(doTimeSeries) {
  for(model in models) {
    if(i %% 100 == 0) {
      print(paste("model:",i))
    }
    i <- i+1
    racipeTS <- subRacipe[,model]
    
    
    
    racipeTS <- sracipeSimulate(racipeTS, numModels = 1,
                                plots = FALSE, genIC = FALSE, anneal = FALSE,
                                genParams = FALSE, integrate = TRUE,
                                integrateStepSize = 0.01,
                                simulationTime = simTime, printInterval = printInterval, printStart = printStart, 
                                initialNoise = iNoise, nNoise = 1, simDet = F, timeSeries = T, stepper = "EMSig")
    # tsMatrix <- t(racipeTS@metadata$timeSeries)
    # tsData_model <- data.frame(tsMatrix, Time=row.names(tsMatrix), Model=model)
    # tsData <- rbind(tsData, tsData_model)
    
    tsData_model <- processDataTS(racipeTS,model, pca, tmpMeans, tmpSds)
    tsData <- rbind(tsData, tsData_model)
    
    
  }
  #sink()
  tsData$Time <- as.numeric(tsData$Time)
  saveRDS(tsData, file = file.path(outputDir,paste0("tsData_",iNoise,"noise_",fc,"FC.Rds")))
  saveRDS(subRacipe, file = file.path(outputDir,paste0("racipeTS_",iNoise,"noise_",fc,"FC.Rds")))
  trialName <- paste0("CTS_transition_",iNoise,"noise_",fc,"FC")
}


####### Plot transitions ####### 
tsData <- readRDS(file.path(outputDir,paste0("tsData_",iNoise,"noise_",fc,"FC.Rds")))

library(gganimate)
library(gifski)
#library(av)
ggplot(tsData, aes(x=-PC1,y=PC2)) +
  geom_point() +
  xlim(-4,4) + ylim(-3,3) +
  labs(title = 'Time: {frame_time}', x = 'PC1', y = 'PC2') +
  transition_time(Time) +
  ease_aes('linear')

anim_save(file.path(getwd(),topoName,"gifs","CTS_clust2_to1_GDx50.gif"))




ggplot() +
  geom_point(data=pca$x, aes(x=PC1,y=PC2)) +
  geom_point(data=tsData[which(tsData$Time == 0.05),], aes(x=PC1,y=PC2),color="red")


ggplot() +
  geom_point(data=pca$x, aes(x=PC1,y=PC2)) +
  geom_point(data=tsData[which(tsData$Time == 60.05),], aes(x=PC1,y=PC2),color="red")

ggplot() +
  geom_point(data=pca$x, aes(x=PC1,y=PC2)) +
  geom_point(data=tsData[which(tsData$Time == 80.00),], aes(x=PC1,y=PC2),color="red")




####### STICCC on snapshot data ####### 
## Load TS data in
racipe_ts <- readRDS(file.path(outputDir,paste0("racipeTS_",iNoise,"noise_",fc,"FC.Rds")))
tsData$Time <- round(tsData$Time, 2)
tsModels <- unique(tsData$Model)
numModels <- length(tsModels)

## Subset data (either downsample models, reduce time granularity, or select one timepoint)
tpInterval <- 0.2
timePoints <- seq(printStart+tpInterval, simTime, tpInterval)
timePoints <- round(timePoints, 2)
#timePoints <- c(5)
dsNum <- numModels


## VICCC in transition
if(useIC) {
  
  processDataIC <- function(racipe,models,pca,tmpMeans,tmpSds) {
    nCols <- nrow(racipe)
    
    normData <- as.data.frame(t(sracipeIC(racipe[,models])))
    normData[,1:nCols] <- log2(1+normData[,1:nCols]) # Log transform
    normData[,1:nCols] <- sweep(normData[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    normData[,1:nCols] <- sweep(normData[,1:nCols], 2, tmpSds, FUN = "/") # scale
    newPca <- (scale(normData[,1:nCols], pca$center, pca$scale) %*%
                 pca$rotation)
    normData$Time <- 0.00
    normData$Model <- models
    normData$PC1 <- newPca[,1]
    normData$PC2 <- newPca[,2]
    normData$PC3 <- newPca[,3]
    
    return(normData)
  }
  
  # Since we only used 1k models, just transfer ICs for them
  ics <- processDataIC(racipe = racipe_ts,
                       models = models,
                       pca = pca,
                       tmpMeans = tmpMeans,
                       tmpSds = tmpSds)
  tsData <- rbind(tsData, ics)
  
} 


## create new stic object for each tpset and calculate velocities
#timepointSets <- list("Untreated"=c(0),"10"=c(10),"20"=c(20),"30"=c(30),"40"=c(40),"50"=c(50),"60"=c(60),"70"=c(70),"80"=c(80))
timepointSets <-  list("Untreated"=c(0),"1"=c(1),"3"=c(3),"5"=c(5),"10"=c(10), "20"=c(20),"40"=c(40),"60"=c(60),"70"=c(70),"80"=c(80))
#timepointSets <-  list("45"=c(45),"50"=c(50),"55"=c(55))
for(tpSet in seq_along(timepointSets)) {
  tpName <- names(timepointSets)[tpSet]
  tps <- timepointSets[tpSet]
  
  metadata <- tsData[which(tsData$Time %in% tps),]
  if(tpName == "80") {
    metadata <- metadata[which(!duplicated(metadata$Model)),]
  }
  metadata$SampleID <- paste0(metadata$Model,"_",metadata$Time)
  rownames(metadata) <- metadata$SampleID
  exprMat <- t(metadata[,genelist])
  
  
  
  sticTS <- sticSE(topo = topo, normData = exprMat,
                topoName = topoName, expName = paste0(topoName, "_transition_t=",tpName,"_2024"))
  
  # sticTS <- SingleCellExperiment(assays = SimpleList(normcounts=exprMat))
  # sticTS@metadata$experimentName <- paste0(topoName,"_time=",tpName,"_Feb22")
  # sticTS@metadata$topoName <- topoName
  # sticTS@metadata$topo <- topo
  # sticTS@metadata$params <- list(sample_radius=0.2, plotScalingFactor=1, gridPlotScalingFactor=1, minNeighbors=5, verbose=T)
  # add metadata
  colData(sticTS) <- DataFrame(metadata)
  colnames(sticTS) <- colData(sticTS)$SampleID
  
  
  # add PCA to SCE object
  posMat <- as.data.frame(scale(t(exprMat), pca$center, pca$scale) %*% pca$rotation)
  rownames(posMat) <- metadata$SampleID
  reducedDim(sticTS, "PCA") <- posMat
  
  sticTS@metadata$pca_data <- pca[1:4]
  sticTS@metadata$pca_summary <- summary(pca)
  
  
  # get grid from previous stic object
  sticTS@metadata$grid.df <- stic@metadata$grid.df
  sticTS@metadata$params$grid.dist <- stic@metadata$params$grid.dist
  sticTS@metadata$params$grid.x.dist <- stic@metadata$params$grid.x.dist
  sticTS@metadata$params$grid.y.dist <- stic@metadata$params$grid.y.dist
  
  # compute pairwise distance between points
  sticTS <- computeDist(sticTS)
  
  # compute trajectories
  sticTS_fname <- file.path(outputDir, paste0("stic_",topoName,"_time=",tpName,"_clust2_GDx50.Rds"))
  
  if(!file.exists(sticTS_fname) | forceSTICCC) {
    sticTS <- runSTICCC(sticTS,  v2=T, invertV2=T)
    saveRDS(sticTS, sticTS_fname)
  } else {
    #sticTS <- readRDS(sticTS_fname) # don't need to read back in for this loop :P
  }
  
}



####### Save/plot results ####### 
## Plot specific timepoints
timepointSets
for(tpSet in seq_along(timepointSets)) {
  tpName <- names(timepointSets)[tpSet]
  
  stic <- readRDS(file.path(outputDir,paste0("stic_",topoName,"_time=",tpName,"_clust2_GDx50.Rds")))
  
  stic <- computeGrid(stic, grid.length = 15)
  stic@metadata$params$xMin <- -3 #stic@metadata$params$xMin + 1
  stic@metadata$params$xMax <- 3 #stic@metadata$params$xMax - 1
  stic@metadata$params$yMin <- -3 #stic@metadata$params$yMin + 1
  stic@metadata$params$yMax <- 3 #stic@metadata$params$yMax - 1
  
  
  
  ## Net
  stic <- computeGridVectors(stic, combine = T, how = "net", unitVectors = F)

  # Flip on the x-axis (PC1) - maybe not necessary?
  colData(stic)$PC1 = -1 * colData(stic)$PC1
  colData(stic)$X = -1 * colData(stic)$X
  stic@metadata$grid.df$dx = -1 * stic@metadata$grid.df$dx
  stic@metadata$grid.df$x.points = -1 * stic@metadata$grid.df$x.points
  
  
  minMagnitude <- 0.001
  scalingFactor <- 1
  arrowheadSize <- 0.3

  plotGrid(sce = stic,
           colorVar = NA,
           plotLoadings = F,
           plotSuffix = paste0("_jul24_net_grey_time=",tpName),
           minMagnitude = minMagnitude,
           scalingFactor = scalingFactor,
           arrowheadSize = arrowheadSize
  )
  
  
  
  ## Rev
  stic <- readRDS(file.path(outputDir,paste0("stic_",topoName,"_time=",tpName,"_clust2_GDx50.Rds")))
  stic <- computeGrid(stic, grid.length = 15)
  stic@metadata$params$xMin <- -3 #stic@metadata$params$xMin + 1
  stic@metadata$params$xMax <- 3 #stic@metadata$params$xMax - 1
  stic@metadata$params$yMin <- -3 #stic@metadata$params$yMin + 1
  stic@metadata$params$yMax <- 3 #stic@metadata$params$yMax - 1
  stic <- computeGridVectors(stic, combine = T, how = "rev", unitVectors = F)
  
  # Flip on the x-axis (PC1) - maybe not necessary?
  colData(stic)$PC1 = -1 * colData(stic)$PC1
  colData(stic)$X = -1 * colData(stic)$X
  stic@metadata$grid.df$dx = -1 * stic@metadata$grid.df$dx
  stic@metadata$grid.df$x.points = -1 * stic@metadata$grid.df$x.points
  
  
  plotGrid(sce = stic,
           colorVar = NA,
           plotLoadings = F,
           plotSuffix = paste0("_jul24_rev_grey_time=",tpName),
           minMagnitude = minMagnitude,
           scalingFactor = scalingFactor,
           arrowheadSize = arrowheadSize
  )
  
  
  
}




####### Basin vs signal strength ####### 

basin_pos <- list()


numTPsPerPhase <- simTime / printInterval / 4
fcs_all <- c(seq(1,50, 50/numTPsPerPhase),rep(50,numTPsPerPhase), seq(50,1, -50/numTPsPerPhase), rep(1,numTPsPerPhase+15))
tps_all <- unique(sort(tsData$Time))

ds_factor <- 50

fcs <- fcs_all[seq(1, length(fcs_all), ds_factor)]
tps <- tps_all[seq(1, length(tps_all), ds_factor)]


for(tp in tps) {
  # Identify basin positions
  pca_df <- tsData[which(tsData$Time == tp),c("PC1","PC2")]
  kde_result <- kde2d(pca_df$PC1, pca_df$PC2, n = 100, lims=c(range(tsData$PC1), range(tsData$PC2)))
  kde_density <- kde_result$z
  energy_landscape <- -log(kde_density)
  
  # Identify local minima as basins
  density_threshold <- 0.1
  local_minima <- matrix(FALSE, nrow = nrow(energy_landscape), ncol = ncol(energy_landscape))
  for (i in 2:(nrow(energy_landscape) - 1)) {
    for (j in 2:(ncol(energy_landscape) - 1)) {
      if (kde_density[i, j] > density_threshold &&
          energy_landscape[i, j] == min(energy_landscape[(i-1):(i+1), (j-1):(j+1)])) {
        local_minima[i, j] <- TRUE
      }
    }
  }
  
  # Convert minima to PCA coordinates
  minima_coords <- which(local_minima, arr.ind = TRUE)
  minima_x <- kde_result$x[minima_coords[, 1]]
  minima_y <- kde_result$y[minima_coords[, 2]]
  
  basin_pos[[as.character(tp)]] <- data.frame(PC1=minima_x, PC2=minima_y, Time=tp)
}

basin_pos_all <- do.call(rbind, basin_pos)


ggplot() +
  geom_point(data=pca$x, aes(x=-PC1,y=PC2),color="grey") +
  geom_point(data=basin_pos_all[,], aes(x=-PC1,y=PC2, color=Time), alpha=0.9) +
  theme_minimal() +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=22),
        axis.line = element_line(color="black", linewidth=0.5),
        axis.ticks = element_line(color="black", linewidth=0.5))



# Track movement of the original minimum
distances <- numeric(length(tps) - 1)
original_min <- basin_pos[["0"]][1, ]  # Take the first minimum from tp=0

for (i in 2:length(tps)) {
  tp <- as.character(tps[i])
  current_minima <- basin_pos[[tp]]
  
  # Compute distances to original minimum
  dists <- sqrt((current_minima$PC1 - original_min$PC1)^2 + (current_minima$PC1 - original_min$PC1)^2)
  
  # Store the minimum distance
  distances[i - 1] <- min(dists)
}

# Store results in a dataframe
distance_df <- data.frame(timepoint = tps[-1], distance = distances)
distance_df <- rbind(data.frame(timepoint=0, distance=0), distance_df)
distance_df$SignalFC <- fcs

ggplot() +
  geom_line(data=distance_df, aes(x=timepoint, y=distance)) +
  xlab("Time") +
  ylab("Distance from original basin") +
  theme_minimal() +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=22))


# Compute scale factor based on data ranges
scale_factor <- max(distance_df$distance, na.rm = TRUE) / max(distance_df$SignalFC, na.rm = TRUE)
ggplot() +
  # First y-axis: Distance
  geom_line(data = distance_df, aes(x = timepoint, y = distance), color = "blue") +
  # Second y-axis: fcs (scaled)
  geom_line(data = distance_df, aes(x = timepoint, y = SignalFC * scale_factor), color = "red") +
  xlab("Time") +
  ylab("Distance from original basin") +
  # Second y-axis label
  scale_y_continuous(
    sec.axis = sec_axis(~ . / scale_factor, name = "SignalFC")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22)
  )


ggplot(distance_df) +
  geom_point(aes(x=SignalFC, y=distance))





