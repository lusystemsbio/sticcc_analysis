##### SETUP AND GLOBALS #####
rm(list=ls())
library(sRACIPE)
library(tidyr)
library(ggplot2)
library(patchwork)
library(STICCC)


# Global parameters
topoName <- "repressilator"
topo <- loadTopo(topoName)
genes <- unique(c(topo$Source, topo$Target))
nGenes <- length(genes)
nEdges <- nrow(topo)
dataDir <- file.path(getwd(),"repressilator","data")
numClusters <- 6
forceSTICCC <- F


## Simulation Parameters
noise <- 0.1
simTime <- 10000
simTimeDet <- simTime/10
printTime <- 0.1
stepSize <- 0.02

# initial conditions
ics <- rep(50, nGenes)


## Model Parameters
g <- 30
k <- 0.5
n <- 4
fc <- 50
th <- 20

##### UTILITY FUNCTIONS #####
source("R/sticcc_analysis_utilities.R")


##### SIMULATION #####
# Set up placeholder racipe object
racipe_sub <- sracipeSimulate(topo, numModels = 1, plots = FALSE, genIC = TRUE,
                              genParams = TRUE, integrate = TRUE, integrateStepSize = stepSize,
                              simulationTime = 0, initialNoise = noise, nNoise = 1, simDet = F, anneal = F)

# Replace parameters with constructed ones
#sracipeParams(racipe_sub)[1,]
sracipeParams(racipe_sub)[1,] <- c(rep(g,nGenes), 
                                   rep(k,nGenes),
                                   rep(th, nEdges), 
                                   rep(n, nEdges), 
                                   rep(fc, nEdges))

# Symmetric ICs
sracipeIC(racipe_sub)[,] <- ics

# Placeholder for zero-noise simulations
racipe_0noise <- racipe_sub

# Long noisy simulation
racipe_fname <- file.path(dataDir, paste0("repressilator_racipe_noise=",noise,"_time=",simTime,"@",printTime))
if(!file.exists(racipe_fname)) {
  racipe_sub <- sracipeSimulate(racipe_sub, simulationTime = simTime,
                                genIC = F, genParams = F,
                                integrateStepSize = stepSize,
                                nNoise = 1, initialNoise = noise,
                                printStart = 0, printInterval = printTime, timeSeries = T,
                                simDet = T, scaledNoise = T)
  saveRDS(racipe_sub, racipe_fname)
} else {
  racipe_sub <- readRDS(racipe_fname)
}



# Zero noise limit cycle estimation
racipe_0noise_fname <- file.path(dataDir, paste0("repressilator_racipe_0noise_time=",simTimeDet,"@",printTime))
if(!file.exists(racipe_0noise_fname)) {
  racipe_0noise <- racipe_sub
  sracipeIC(racipe_0noise)[,] <- c(50, 40, 45)
  racipe_0noise <- sracipeSimulate(racipe_0noise, simulationTime = simTimeDet,
                                   genIC = F, genParams = F,
                                   integrateStepSize = stepSize,
                                   nNoise = 1, initialNoise = 0.1,
                                   printStart = 0, printInterval = printTime, timeSeries = T,
                                   simDet = T, scaledNoise = T)
  saveRDS(racipe_0noise, racipe_0noise_fname)
} else {
  racipe_0noise <- readRDS(racipe_0noise_fname)
}




##### TRAJECTORY PLOTTING #####
# Extract trajectory (wide)
traj_dat <- as.data.frame(t(racipe_sub@metadata$timeSeries))
colnames(traj_dat) <- genes
traj_dat$Time <- seq(printTime,simTime+printTime,printTime)
for(gene in genes) {
  traj_dat[,gene] <- as.numeric(traj_dat[,gene])  
}


# Plot expression of all genes over time (needs trajectory in long format)
traj_df <- pivot_longer(traj_dat, cols = all_of(genes), names_to = "Gene", values_to = "Expression")
traj_df$Expression <- as.numeric(traj_df$Expression)

traj_plot <- ggplot(traj_df[which(traj_df$Time < 100),], aes(x=Time,y=Expression, color=as.factor(Gene))) +
  geom_point() 
traj_plot


##### PLOTTING ON WT PCA #####
# Plot trajectory on WT racipe PCA
topoName <- "repressilator"
topoDir <- file.path(getwd(),topoName)
outputDir = file.path(topoDir,"data")
plotDir <- file.path(topoDir,"gridpt_val")
if(!dir.exists(plotDir)) {
  dir.create(plotDir, recursive = T)
}


racipe.ensemble <- readRDS(file.path(outputDir, paste0("simData_",topoName,".Rds")))
pca <- readRDS(file.path(outputDir,paste0("PCA_res.Rds")))
simExp <- assay(racipe.ensemble, 1)
simExp <- log2(1+simExp)
tmpMeans <- rowMeans(simExp)
tmpSds <- apply(simExp,1,sd)

# normalize trajectory data
nCols <- length(genes)
traj_norm <- as.data.frame(traj_dat[,c(1:nCols)])
traj_norm[,1:nCols] <- log2(1+traj_norm[,1:nCols]) # Log transform
traj_norm[,1:nCols] <- sweep(traj_norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
traj_norm[,1:nCols] <- sweep(traj_norm[,1:nCols], 2, tmpSds, FUN = "/") # scale

traj_pca <- as.data.frame(predict(pca, traj_norm))
traj_pca$Time <- traj_dat$Time


loadingDF <- as.data.frame(pca$rotation)
loadingDF$gene <- rownames(loadingDF)
loadingDF$X <- 0
loadingDF$Y <- 0


loadingFactor <- 3
loadingLabelFactor <- loadingFactor + 0.5
ggplot(data=traj_pca[which(traj_pca$Time > 9979),], aes(x=PC1,y=PC2, color=as.numeric(Time))) + 
  geom_point(data=as.data.frame(pca$x), aes(x=PC1,y=PC2), color="grey", alpha=0.5) +
  geom_point(aes(color=as.numeric(Time))) +
  geom_path(arrow = arrow(length = unit(0.075, "inches"))) + 
  geom_segment(data = loadingDF[,], aes(xend=PC1*loadingFactor, yend=PC2*loadingFactor, x = X, y = Y), color="black", arrow = arrow(length = unit(0.6,"cm"))) +
  geom_text(data = loadingDF[,], aes(x = PC1*loadingLabelFactor, y = PC2*loadingLabelFactor, label=gene), color="black", size=8, fontface="bold") +
  ggtitle(paste0("Noise=",noise,"_time=",simTime,"@",printTime))



##### STICCC #####
# Create sticcc object, partition data into grid points
# Subset 2k points to use for STICCC
subset_fname <- file.path(dataDir,"validation_subset_ids_07302024.Rds")
if(!file.exists(subset_fname)) {
  vector_subset <- sample(1:nrow(traj_norm), 2000)
  saveRDS(vector_subset, subset_fname)
} else {
  vector_subset <- readRDS(subset_fname)
}


exprMat_norm <- as.matrix(t(traj_norm[vector_subset,]))
stic <- sticSE(topo = topo, normData = exprMat_norm,
               topoName = topoName, expName = paste0("grid_val_noise=",noise,"_time=",simTime,"@",printTime))


# add metadata
stic <- prepMetadata(stic, exprMat_norm, cluster = T, k = numClusters)
stic@metadata$params$sample_radius <- 0.15
colData(stic)$Time <- traj_pca[vector_subset,"Time"]


# add PCA 
reducedDim(stic, "PCA") <- traj_pca[vector_subset,-which(colnames(traj_pca) == "Time")]
stic@metadata$params$nPCs <- ncol(traj_pca)-1


# set borders if backdrop provided
borders <- list(xmin=floor(min(traj_pca[,1])),xmax=ceiling(max(traj_pca[,1])),ymin=floor(min(traj_pca[,2])),ymax=ceiling(max(traj_pca[,2])))
stic@metadata$params$xMin <- borders[['xmin']]
stic@metadata$params$xMax <- borders[['xmax']]
stic@metadata$params$yMin <- borders[['ymin']]
stic@metadata$params$yMax <- borders[['ymax']]


# compute grid based on PCA
stic <- computeGrid(stic, xMin = -3, xMax=3, yMin=-3, yMax=3)

# compute pairwise distance between points
stic <- computeDist(stic)

## Compute vectors
stic_fname <- file.path(outputDir, paste0("stic_",stic@metadata$experimentName,"_subset2k.Rds"))
if(!file.exists(stic_fname) | forceSTICCC) {
  stic <- runSTICCC(stic, v2=T, invertV2=T)
  saveRDS(stic, stic_fname)
} else {
  stic <- readRDS(stic_fname)
}


##### LIMIT CYCLE VALIDATION #####
# Extract deterministic limit cycle & annotate time
traj_0noise <- as.data.frame(t((racipe_0noise@metadata$timeSeriesDet)))
traj_0noise$Time <- round(seq(printTime,simTimeDet+printTime,printTime), 1)
for(gene in genes) {
  traj_0noise[,gene] <- as.numeric(traj_0noise[,gene])  
}


# Normalize with moments from original simulation
traj_0noise_norm <- traj_0noise[,c(1:nCols)]
traj_0noise_norm[,1:nCols] <- log2(1+traj_0noise_norm[,1:nCols]) # Log transform
traj_0noise_norm[,1:nCols] <- sweep(traj_0noise_norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
traj_0noise_norm[,1:nCols] <- sweep(traj_0noise_norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
colnames(traj_0noise_norm) <- genes

# Project to same PCs
traj_0noise_pca <- as.data.frame(predict(pca, traj_0noise_norm))
traj_0noise_pca$Time <- round(traj_0noise$Time,1)



# Get data for one limit cycle
cycleSize <- 7.6
cycleStart <- 994.7-cycleSize*10
cycleEnd <- cycleStart + cycleSize
limitCycleTimes <- seq(cycleStart, cycleEnd, printTime)
limitCycleData <- traj_0noise[which(traj_0noise$Time %in% limitCycleTimes),
                              which(!colnames(traj_0noise_pca) == "Time")]


# Plot one limit cycle (sanity check)
ggplot(traj_0noise_pca[which(traj_0noise_pca$Time > cycleStart & traj_0noise_pca$Time < cycleEnd),], aes(x=PC1,y=PC2,color=Time)) +
  geom_point() +
  geom_path()


# Compute vectors for each point along a limit cycle
#cycleNo <- 0
#cycleEnd <- 1000.1 - cycleNo*cycleSize
#cycleStart <- cycleEnd - cycleSize

cycleTs <- round(seq(cycleStart, cycleEnd, printTime), 1)
cycleData <- traj_0noise_norm[which(traj_0noise_pca$Time %in% round(cycleTs,1)),
                              which(!colnames(traj_0noise_pca) == "Time")]
rownames(cycleData) <- cycleTs

cycleDataPCA <- traj_0noise_pca[which(traj_0noise_pca$Time %in% round(cycleTs,1)),
                                which(!colnames(traj_0noise_pca) == "Time")]
rownames(cycleDataPCA) <- cycleTs

v_df <- data.frame(Time=cycleTs,
                   v1x=NA,v1y=NA,
                   v2x=NA,v2y=NA,
                   netx=NA,nety=NA,
                   revx=NA,revy=NA)


for(t in seq_along(cycleTs)) {
  v_smooth <- smoothVector(sce = stic,
                           queryPoint = cycleDataPCA[t,],
                           neighborhoodRadius = 0.1,
                           invertV2 = T)
  if(!is.null(v_smooth)) {
    v_df$v1x[t] <- v_smooth$v1_x
    v_df$v1y[t] <- v_smooth$v1_y
    
    v_df$v2x[t] <- v_smooth$v2_x
    v_df$v2y[t] <- v_smooth$v2_y
    
    v_df$netx[t] <- v_smooth$net_x
    v_df$nety[t] <- v_smooth$net_y
    
    v_df$revx[t] <- v_smooth$rev_x
    v_df$revy[t] <- v_smooth$rev_y
  } else {
    v_df$v1x[t] <- 0
    v_df$v1y[t] <- 0
    
    v_df$v2x[t] <- 0
    v_df$v2y[t] <- 0
    
    v_df$netx[t] <- 0
    v_df$nety[t] <- 0
    
    v_df$revx[t] <- 0
    v_df$revy[t] <- 0
  }
}

scale_factor <- 0.5
v_df$x <- traj_0noise_pca[which(traj_0noise_pca$Time %in% cycleTs), "PC1"]
v_df$y <- traj_0noise_pca[which(traj_0noise_pca$Time %in% cycleTs), "PC2"]
ggplot(v_df, aes(x=x,y=y,color=Time)) +
  geom_point() +
  geom_segment(aes(xend=x+v1x*scale_factor, yend=y+v1y*scale_factor), color="black", arrow = arrow(length = unit(0.4,"cm")))


image <- ggplot() + 
  geom_point(data=as.data.frame(pca$x), aes(x=PC1,y=PC2), color="grey", alpha=0.5) +
  geom_point(data=traj_pca[which(traj_pca$Time > 9959),], aes(x=PC1,y=PC2, color=Time)) +
  geom_path(data=traj_pca[which(traj_pca$Time > 9959),], aes(x=PC1,y=PC2,color=Time), arrow = arrow(length = unit(0.075, "inches"))) + 
  geom_segment(data = loadingDF[,], aes(xend=PC1*loadingFactor, yend=PC2*loadingFactor, x = X, y = Y), color="black", arrow = arrow(length = unit(0.6,"cm"))) +
  geom_text(data = loadingDF[,], aes(x = PC1*loadingLabelFactor, y = PC2*loadingLabelFactor, label=gene), color="black", size=8, fontface="bold") +
  ggtitle(paste0("Noise=",noise,"_time=",simTime,"@",printTime)) +
  geom_point(data=v_df, aes(x=x,y=y), color="black") #+
  #geom_segment(data=v_df, aes(x=x, y=y, xend=x+v1x*scale_factor, yend=y+v1y*scale_factor), color="black", arrow = arrow(length = unit(0.4,"cm")))
image

fname <- file.path(plotDir,paste0("PCA_noisy_and_det_trajectories_loadings.pdf"))
pdf(fname, width = 10, height = 10)
print(image)
dev.off()


####### VALIDATION W/ OPTIMIZED LAG #######
# Take a subset of points along the limit cycle
subset_trajectory <- cycleDataPCA[seq(2,74,8),]
rownames(subset_trajectory) <- seq(2,74,8)

# Calculate vectors for each subsetted point
traj_v_pred <- trajectorySmoothVectors(trajectory = subset_trajectory, # trajectory in PCA coordinates
                                       sce = stic,
                                       neighborhoodRadius = 0.05,
                                       invertV2 = T,
                                       vec.use = "net") 


traj_v_pred_REV <- trajectorySmoothVectors(trajectory = subset_trajectory, # trajectory in PCA coordinates
                                           sce = stic,
                                           neighborhoodRadius = 0.05,
                                           invertV2 = T,
                                           vec.use = "rev") 

## Identify lag for each point to keep variance constant
var_traj_df <- data.frame()
#optimize_lags <- seq(0.1, 3, 1)
optimize_lags <- c(seq(0.1, 2.9, 0.1), seq(3, 20, 1))

for(subsetPtIdx in rownames(subset_trajectory)) {
  subsetPt <- subset_trajectory[subsetPtIdx,]
  
  ptVar <- vobs_var_by_t(trajectory = traj_pca,
                         lags = optimize_lags,
                         sce = stic,
                         queryPoint = subsetPt,
                         neighborhoodRadius = 0.02,
                         plot=F, 
                         save=F)
  
  ptVar$QueryPoint <- subsetPtIdx
  var_traj_df <- rbind(var_traj_df, ptVar)
  
}

var_traj_df$QueryPoint <- factor(var_traj_df$QueryPoint, levels=as.character(sort(as.numeric(rownames(subset_trajectory)))))

# Look at variance over time
ggplot(var_traj_df) +
  geom_point(aes(x=Lag, y=Var, color=QueryPoint), size=3) +
  geom_line(aes(x=Lag, y=Var, color=QueryPoint)) +
  theme_sticcc()

# RMSD between initial and final over time
ggplot(var_traj_df) +
  geom_point(aes(x=Lag, y=RMSD, color=QueryPoint), size=3) +
  geom_line(aes(x=Lag, y=RMSD, color=QueryPoint)) +
  theme_sticcc()






# identify lags closest to a specified threshold RMSD
#target_rmsd <- 0.2
target_rmsd <- max(var_traj_df[which(var_traj_df$Lag == 0.1), "RMSD"])

optimal_lags <- data.frame(QueryPoint = unique(var_traj_df$QueryPoint), OptimalLag = NA, Var = NA)
for(pt in unique(var_traj_df$QueryPoint)) {
  # Find lag where var is closest to target_var
  pt_vars <- var_traj_df[which(var_traj_df$QueryPoint == pt & var_traj_df$Lag <= 10),]
  pt_lag <- pt_vars[which.min(abs(pt_vars$RMSD - target_rmsd)), "Lag"]
  pt_minVar <- pt_vars[which.min(abs(pt_vars$RMSD - target_rmsd)), "RMSD"]
  
  optimal_lags[which(optimal_lags$QueryPoint == pt), "OptimalLag"] <- pt_lag
  optimal_lags[which(optimal_lags$QueryPoint == pt), "RMSD"] <- pt_minVar
}

#plot(optimal_lags$QueryPoint, optimal_lags$OptimalLag)


ggplot(var_traj_df[which(var_traj_df$Lag < 2),]) +
  geom_point(aes(x=Lag, y=RMSD, color=QueryPoint), size=3) +
  geom_line(aes(x=Lag, y=RMSD, color=QueryPoint)) +
  theme_sticcc() +
  geom_hline(aes(yintercept=target_rmsd), color="red", lty="dashed", size=2) +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))




#debug(v_obs_along_path)
rs_list <- v_obs_along_path(trajectory = traj_pca, # PCA plus Time column
                            lag = optimal_lags$OptimalLag, # numeric - time gap (computed numerically, so should correspond to times, not indices!)
                            sce = stic,
                            queryTrajectory = subset_trajectory, # either rownames of trajectory, or a dataframe of same ncol to be compared to it
                            neighborhoodRadius = 0.02,
                            v_pred = traj_v_pred) 



rs_summary <- rs_list$Summary
rs_boxplot <- rs_list$BoxplotData
rs_boxplot$QueryPoint <- factor(rs_boxplot$QueryPoint, levels=as.character(seq(2,74,8)))



# Scale angles and generate vectors
vObs_scaling_factor <- 0.5
rs_summary$Obs.Mag.1.Scaled <- rs_summary$Obs.Mag.1 / max(rs_summary$Obs.Mag.1) * vObs_scaling_factor
rs_summary$Obs.Mag.2.Scaled <- rs_summary$Obs.Mag.2 / max(rs_summary$Obs.Mag.2) * vObs_scaling_factor

rs_summary$Obs.Vector.X.1 <- cos(rs_summary$Obs.Angle.1) * rs_summary$Obs.Mag.1.Scaled
rs_summary$Obs.Vector.Y.1 <- sin(rs_summary$Obs.Angle.1) * rs_summary$Obs.Mag.1.Scaled

rs_summary$Obs.Vector.X.2 <- cos(rs_summary$Obs.Angle.2) * rs_summary$Obs.Mag.2.Scaled
rs_summary$Obs.Vector.Y.2 <- sin(rs_summary$Obs.Angle.2) * rs_summary$Obs.Mag.2.Scaled


# Using sampled ideal path, calculate a fwd and back angle for each point (exclude first and last)
hline_df <- data.frame(QueryPoint=rs_summary$QueryPoint,
                       Angle.Fwd = NA,
                       Angle.Back = NA,
                       Angle.NetFlow = NA,
                       Angle.Rev = NA
)
for(i in 1:nrow(subset_trajectory)) {
  angle_fwd <- NA
  angle_back <- NA
  
  # Calculate angle to previous row
  if(i > 1) {
    diff_back <- subset_trajectory[i-1,c("PC1","PC2")] - subset_trajectory[i,c("PC1","PC2")]
    angle_back <- angle_conversion(diff_back)
    
  }
  # Angle to next row
  if(i < nrow(subset_trajectory)) {
    diff_fwd <- subset_trajectory[i+1,c("PC1","PC2")] - subset_trajectory[i,c("PC1","PC2")]
    angle_fwd <- angle_conversion(diff_fwd)
  }
  
  # Store output
  hline_df[i,"Angle.Fwd"] <- angle_fwd
  hline_df[i,"Angle.Back"] <- angle_back
  hline_df[i,"Angle.NetFlow"] <- angle_conversion(traj_v_pred[i, c("dx","dy")]) 
  hline_df[i,"Angle.Rev"] <- angle_conversion(traj_v_pred_REV[i, c("dx","dy")]) 
  hline_df[i,"Angle.Rev.180"] <- angle_conversion(-1*traj_v_pred_REV[i, c("dx","dy")]) 
  
  
}


blue_segment_df <- data.frame(x=1, xend=9, y=hline_df$Angle.Fwd[1], yend=hline_df$Angle.Fwd[9])


# Violin plot with annotations
image <- ggplot(data=rs_boxplot) +
  geom_violin(aes(x=QueryPoint, y=Angle)) +
  #geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.Fwd), color = "blue", size = 4, alpha=0.8) +
  geom_segment(data = blue_segment_df, aes(x = x, xend = xend, y = y, yend=yend), color = "blue", size = 0.7, linetype="dashed", alpha=0.8) +
  geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.NetFlow), color = "purple", size = 4, alpha=0.8) +
  labs(x = "Trajectory Point") +
  theme_sticcc()
image

pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD.pdf")), width = 10, height = 10)
print(image)
dev.off()


#library(tidyverse)
library(viridisLite)

mywidth <- .45 # bit of trial and error
p <- ggplot(rs_boxplot) + geom_violin(aes(x=QueryPoint,y=Angle))

# all you need for the gradient fill
vl_fill <- data.frame(ggplot_build(p)$data) %>%
  mutate(xnew = x - mywidth * violinwidth, xend = x + mywidth * violinwidth)

breaks <- unique(as.integer(rs_boxplot$QueryPoint))
labels <- unique(rs_boxplot$QueryPoint)

image <- ggplot() +
  geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y,
                                   color = violinwidth), show.legend = FALSE) +
  # Re-use geom_violin to plot the outline
  geom_violin(data = rs_boxplot, aes(x = as.integer(QueryPoint), y = Angle, fill = QueryPoint),
              color = "white", alpha = 0, draw_quantiles = c(0.25, 0.5, 0.75),
              show.legend = FALSE) +
  scale_x_continuous(breaks = breaks, labels = labels) +
  scale_color_viridis_c() +
  geom_segment(data = blue_segment_df, aes(x = x, xend = xend, y = y, yend=yend), color = "blue", size = 0.7, linetype="dashed", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:10, y = Angle.NetFlow), color = "purple", size = 4, alpha=0.8) +
  theme_sticcc() +
  labs(x = "Trajectory Point", y = "Angle")
image

pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill.pdf")), width = 10, height = 10)
print(image)
dev.off()



# Density plot with reversibility
revScalingFactor <- 0.1
image <- ggplot() +
  #geom_point(data=pca$x, aes(x=PC1, y=PC2)) +
  geom_density2d(data=traj_pca[vector_subset,], aes(x=PC1, y=PC2), color="red") +
  geom_point(data=cycleDataPCA, aes(x=PC1, y=PC2, color=1:nrow(cycleDataPCA))) +
  geom_segment(data = traj_v_pred, 
               aes(x=x,y=y, xend=x+dx*revScalingFactor, yend=y+dy*revScalingFactor), 
               arrow = arrow(length = unit(0.3,"cm")), color="black", size=2, alpha=0.7) +
  scale_color_gradient(name="Point No.", breaks=c(20,40,60)) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))


pdf(file = file.path(plotDir, paste0("VObs_On_PCA.pdf")), width = 10, height = 10)
print(image)
dev.off()






# Example plot of init and final states
selectedPtIdx <- 58
selectedPt <-  cycleDataPCA[selectedPtIdx,]

image <- plot_v_vs_trajectory(sce = stic,
                              queryPoint = selectedPt, 
                              neighborhoodRadius = 0.02,
                              trajectory = traj_pca,
                              lag = 0.5,
                              v_pred = traj_v_pred[8,], 
                              scale_factor = 0.2)
image + theme_sticcc() + theme(axis.line = element_line(linewidth = 0.7, colour = "black"))



image <- ggplot() +
  #geom_point(data=pca$x, aes(x=PC1, y=PC2)) +
  geom_density2d(data=traj_pca[vector_subset,], aes(x=PC1, y=PC2), color="red") +
  geom_point(data=traj_0noise_pca[which(traj_0noise_pca$Time %in% cycleTs),], aes(x=PC1, y=PC2, color=1:77)) +
  geom_segment(data = traj_v_pred, 
               aes(x=x,y=y, xend=x+dx*0.2, yend=y+dy*0.2), 
               arrow = arrow(length = unit(0.3,"cm")), color="black", size=2, alpha=0.7) +
  scale_color_gradient(name="Query Point", breaks=c(10,30,50)) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))
image



## Extended data figure to show deterministic limit cycle and point indices
image <- ggplot() +
  #geom_point(data=pca$x, aes(x=PC1, y=PC2)) +
  geom_density2d(data=traj_pca[vector_subset,], aes(x=PC1, y=PC2), color="blue") +
  geom_point(data=traj_0noise_pca[which(traj_0noise_pca$Time %in% cycleTs),], aes(x=PC1, y=PC2, color=1:77)) +
  #geom_segment(data = traj_v_pred, 
  #             aes(x=x,y=y, xend=x+dx*0.2, yend=y+dy*0.2), 
  #             arrow = arrow(length = unit(0.3,"cm")), color="black", size=2, alpha=0.7) +
  scale_color_gradient(name="Query Point", breaks=c(10,30,50)) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))
image




















