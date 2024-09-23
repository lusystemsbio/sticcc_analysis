##### SETUP AND GLOBALS #####
rm(list=ls())
library(sRACIPE)
library(tidyr)
library(ggplot2)
library(STICCC)

# A |----| B
# |        ^
# V        |
# C |----| D

##### UTILITY FUNCTIONS #####
source("R/sticcc_analysis_utilities.R")

# Feb 2024: consider looking into the behavior of the topology with B-->D swapped in direction (simpler circuit w/o complete feedback, likely)
# Main states should be A/C and B/D
topo <- data.frame(Source=c("A","A","B","C","D","D"), 
                   Target=c("B","C","A","D","B","C"), 
                   Type=c(2,1,2,2,1,2))
genes <- unique(c(topo$Source, topo$Target))
nGenes <- length(genes)
nEdges <- nrow(topo)
dataDir <- file.path(getwd(),"CTS","data")
set.seed(123)
forceSTICCC <- F


## Simulation Parameters
noise <- 2
simTime <- 100000
printTime <- 0.1
stepSize <- 0.02
nParts <- 10
numClusters <- 4

# initial conditions
ics <- rep(50, nGenes)


## Model Parameters
g <- 50
k <- 0.1
n <- 4
fc <- 10
# modified threshold is lower to make B --| A and C --| D more active 
th <- 100
th_mod <- 80




##### SIMULATION #####
# Set up placeholder racipe object
racipe_sub <- sracipeSimulate(topo, numModels = 1, plots = FALSE, genIC = TRUE,
                              genParams = TRUE, integrate = TRUE, integrateStepSize = stepSize,
                              simulationTime = 0, initialNoise = noise, nNoise = 1, simDet = F, anneal = F)

# Replace parameters
sracipeParams(racipe_sub)[1,]
sracipeParams(racipe_sub)[1,] <- c(rep(g,nGenes), 
                                   rep(k,nGenes),
                                   c(th_mod,rep(th, 4),th_mod), # B_A, A_B, D_B, A_C, D_C, C_D
                                   rep(n, nEdges), 
                                   rep(fc, nEdges))

# Symmetric ICs
sracipeIC(racipe_sub)[,] <- ics

# Long noisy simulation
for(part in 1:nParts) {
  racipe_fname <- file.path(dataDir,paste0("CTS_racipe_noise=",noise,"_time=",simTime,"@",printTime,"_part=",part))
  if(!file.exists(racipe_fname)) {
    if(part == 1) {
      racipe_sub <- sracipeSimulate(racipe_sub, simulationTime = simTime/nParts,
                                    genIC = F, genParams = F,
                                    integrateStepSize = stepSize,
                                    nNoise = 1, initialNoise = noise,
                                    printStart = 0, printInterval = printTime, timeSeries = T,
                                    simDet = T, scaledNoise = T)
    } else {
      prev_fname <- file.path(dataDir,paste0("CTS_racipe_noise=",noise,"_time=",simTime,"@",printTime,"_part=",(part-1)))
      prev_part <- readRDS(prev_fname)
      
      # Continue where we left off
      sracipeIC(racipe_sub) <- unlist(prev_part@metadata$timeSeries[,ncol(prev_part@metadata$timeSeries)])
      racipe_sub <- sracipeSimulate(racipe_sub, simulationTime = simTime/nParts,
                                    genIC = F, genParams = F,
                                    integrateStepSize = stepSize,
                                    nNoise = 1, initialNoise = noise,
                                    printStart = 0, printInterval = printTime, timeSeries = T,
                                    simDet = T, scaledNoise = T)
    }
    
    saveRDS(racipe_sub, racipe_fname)
  } 
}


##### JOINING SIMULATION PARTS #####
tsData <- data.frame()
for(part in 1:nParts) {
  racipe_part <- readRDS(file.path(dataDir,paste0("CTS_racipe_noise=",noise,"_time=",simTime,"@",printTime,"_part=",part)))
  tsData_part <- as.data.frame(t(racipe_part@metadata$timeSeries))
  tsData <- rbind(tsData, tsData_part)
}




##### TRAJECTORY PLOTTING #####
# Extract trajectory (wide)
traj_dat <- tsData
colnames(traj_dat) <- genes
traj_dat$Time <- seq(printTime,simTime+printTime*nParts,printTime)
for(gene in genes) {
  traj_dat[,gene] <- as.numeric(traj_dat[,gene])  
}
#traj_dat$State <-   ifelse(((traj_dat$A + traj_dat$C) > (traj_dat$B + traj_dat$D)), "A/C", "B/D")



# Plot expression of all genes over time (needs trajectory in long format)
traj_df <- pivot_longer(traj_dat, cols = all_of(genes), names_to = "Gene", values_to = "Expression")
traj_df$Expression <- as.numeric(traj_df$Expression)

traj_plot <- ggplot(traj_df[which(traj_df$Time < 1000),], aes(x=Time,y=Expression, color=as.factor(Gene))) +
  geom_point() 
traj_plot


##### PLOTTING ON WT PCA #####
# Plot trajectory on WT racipe PCA
topoName <- "CTS"
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
ggplot(data=traj_pca[which(traj_pca$Time > 99979),], aes(x=PC1,y=PC2, color=as.numeric(Time))) + 
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
  vector_subset <- sample(1:nrow(traj_norm), 5000)
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
stic_fname <- file.path(outputDir, paste0("stic_",stic@metadata$experimentName,"_subset5k.Rds"))
if(!file.exists(stic_fname) | forceSTICCC) {
  stic <- runSTICCC(stic, v2=T, invertV2=T)
  saveRDS(stic, stic_fname)
} else {
  stic <- readRDS(stic_fname)
}



##### DEFINE LINEAR PATH #####
## Use linear interpolation to create a set of points from the center of cluster 1 to that of cluster 2, via cluster 3 (or 4)
## Ideal trajectory should be in PCA space
traj_sample_fname <- file.path(dataDir,"validation_trajectory_sample.Rds")
if(!file.exists(traj_sample_fname)) {
  traj_sample <- c(vector_subset, sample(1:nrow(traj_pca), 10000))
  saveRDS(traj_sample, traj_sample_fname)
} else {
  traj_sample <- readRDS(traj_sample_fname)
}

traj_clust_fname <- file.path(dataDir,"validation_trajectory_clusters.Rds")
if(!file.exists(traj_clust_fname)) {
  clusters <- kmeans(traj_pca[traj_sample,c(1:4)], centers=3)$cluster 
  saveRDS(clusters, traj_clust_fname)
} else {
  clusters <- readRDS(traj_clust_fname)
}

#ggplot(pca$x, aes(x=PC1,y=PC2)) +
#  geom_point(aes(color=as.factor(clusters)))

#ggplot(pca$x, aes(x=PC1,y=PC3)) +
#  geom_point(aes(color=as.factor(clusters)))

ggplot(traj_pca[traj_sample,], aes(x=PC1,y=PC2)) +
  geom_point(aes(color=as.factor(clusters)))

clust4_centroid <- colMeans(traj_pca[traj_sample[which(clusters == 1)],colnames(pca$x)])
clust3_centroid <- colMeans(traj_pca[traj_sample[which(clusters == 2)],colnames(pca$x)])
clust1_centroid <- colMeans(traj_pca[traj_sample[which(clusters == 3)],colnames(pca$x)])


ggplot(traj_pca[traj_sample,], aes(x=PC1,y=PC2)) +
  geom_point(aes(color=as.factor(clusters))) +
  geom_point(data=t(as.data.frame(clust4_centroid)), aes(x=PC1, y=PC2), size=3, color="black") +
  geom_point(data=t(as.data.frame(clust3_centroid)), aes(x=PC1, y=PC2), size=3, color="black") +
  geom_point(data=t(as.data.frame(clust1_centroid)), aes(x=PC1, y=PC2), size=3, color="black")

clust_431_traj <- rbind(interpolate_vectors(clust4_centroid, clust3_centroid, 30), 
                        interpolate_vectors(clust3_centroid, clust1_centroid, 30))
colnames(clust_431_traj) <- colnames(pca$x)


ggplot() +
  #geom_point(data=pca$x, aes(x=PC1, y=PC2)) +
  geom_density2d(data=traj_pca[vector_subset,], aes(x=PC1, y=PC2)) +
  geom_point(data=clust_431_traj, aes(x=PC1, y=PC2, color=1:60)) +
  scale_color_gradient(name="Query Point") +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))




##### STICCC PLOTS? #####

# Plot results
minMagnitude <- 0.001
scalingFactor <- 0.5
arrowheadSize <- 0.4

colData(stic)$Cluster <- clusters[1:5000]




##### V_OBS ALONG IDEAL PATH #####
subset_trajectory <- clust_431_traj[c(1,seq(6,60,6)),]


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
optimize_lags <- c(seq(0.1, 2.9, 0.1), seq(3, 80, 1))

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

#plotPt <- rownames(subset_trajectory)[1]
ggplot(var_traj_df) +
  geom_point(aes(x=Lag, y=Var, color=QueryPoint), size=3) +
  geom_line(aes(x=Lag, y=Var, color=QueryPoint)) +
  theme_sticcc()

ggplot(var_traj_df[which(var_traj_df$Lag < 10),]) +
  geom_point(aes(x=Lag, y=RMSD, color=QueryPoint), size=3) +
  geom_line(aes(x=Lag, y=RMSD, color=QueryPoint)) +
  theme_sticcc()




# identify lags closest to a specified threshold variance
# Select target rmsd by the maximum rmsd @ lag=0.01
target_rmsd <- max(var_traj_df[which(var_traj_df$Lag == 0.1), "RMSD"])
optimal_lags <- data.frame(QueryPoint = unique(var_traj_df$QueryPoint), OptimalLag = NA, Var = NA)
for(pt in unique(var_traj_df$QueryPoint)) {
  # Find lag where var is closest to target_var
  pt_vars <- var_traj_df[which(var_traj_df$QueryPoint == pt),]
  pt_lag <- pt_vars[which.min(abs(pt_vars$RMSD - target_rmsd)), "Lag"]
  pt_minVar <- pt_vars[which.min(abs(pt_vars$RMSD - target_rmsd)), "RMSD"]
  
  optimal_lags[which(optimal_lags$QueryPoint == pt), "OptimalLag"] <- pt_lag
  optimal_lags[which(optimal_lags$QueryPoint == pt), "RMSD"] <- pt_minVar
}

#plot(optimal_lags$QueryPoint, optimal_lags$OptimalLag)

ggplot(var_traj_df[which(var_traj_df$Lag < 10),]) +
  geom_point(aes(x=Lag, y=RMSD, color=QueryPoint), size=3) +
  geom_line(aes(x=Lag, y=RMSD, color=QueryPoint)) +
  theme_sticcc() +
  geom_hline(aes(yintercept=target_rmsd), color="red", lty="dashed", size=2) +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))




rs_list <- v_obs_along_path(trajectory = traj_pca, # PCA plus Time column
                            lag = optimal_lags$OptimalLag, # numeric - time gap (computed numerically, so should correspond to times, not indices!)
                            sce = stic,
                            queryTrajectory = subset_trajectory, # either rownames of trajectory, or a dataframe of same ncol to be compared to it
                            neighborhoodRadius = 0.02,
                            v_pred = traj_v_pred) 



rs_summary <- rs_list$Summary
rs_boxplot <- rs_list$BoxplotData
rs_boxplot$QueryPoint <- factor(rs_boxplot$QueryPoint, levels=c("1","6","12","18","24","30","36","42","48","54","60"))


# Scale angles and generate vectors
vObs_scaling_factor <- 0.5
rs_summary$Obs.Mag.1.Scaled <- rs_summary$Obs.Mag.1 / max(rs_summary$Obs.Mag.1) * vObs_scaling_factor
rs_summary$Obs.Mag.2.Scaled <- rs_summary$Obs.Mag.2 / max(rs_summary$Obs.Mag.2) * vObs_scaling_factor

rs_summary$Obs.Vector.X.1 <- cos(rs_summary$Obs.Angle.1) * rs_summary$Obs.Mag.1.Scaled
rs_summary$Obs.Vector.Y.1 <- sin(rs_summary$Obs.Angle.1) * rs_summary$Obs.Mag.1.Scaled

rs_summary$Obs.Vector.X.2 <- cos(rs_summary$Obs.Angle.2) * rs_summary$Obs.Mag.2.Scaled
rs_summary$Obs.Vector.Y.2 <- sin(rs_summary$Obs.Angle.2) * rs_summary$Obs.Mag.2.Scaled





## Violin plot with annotations:
# Black: dir fwd and backward on ideal path
# Purple: net flow
# Purple dashed: reversibility



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




image <- ggplot(data=rs_boxplot) +
  geom_violin(aes(x=QueryPoint, y=Angle)) +
  #geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.Fwd), color = "blue", size = 4, alpha=0.8) +
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("1","6","12","18","24")),], aes(x = 1, xend = 5, y = Angle.Fwd, yend = Angle.Fwd), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("1","6","12","18","24")),], aes(x = 1, xend = 5, y = Angle.Back, yend = Angle.Back), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("36","42","48","54","60")),], aes(x = 7, xend = 11, y = Angle.Fwd, yend = Angle.Fwd), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("36","42","48","54","60")),], aes(x = 7, xend = 11, y = Angle.Back, yend = Angle.Back), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  #geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.Back), color = "blue", size = 4, , alpha=0.8) +
  geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.Rev), color = "purple", size = 4, alpha=0.8) +
  geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.Rev.180), color = "purple", size = 4, alpha=0.8) +
  labs(x= "Trajectory Point") +
  theme_sticcc()
image


pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD.pdf")), width = 10, height = 10)
print(image)
dev.off()

# Now plot with gradient fill on violin plot
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
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("1","6","12","18","24")),], aes(x = 1, xend = 5, y = Angle.Fwd, yend = Angle.Fwd), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("1","6","12","18","24")),], aes(x = 1, xend = 5, y = Angle.Back, yend = Angle.Back), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("36","42","48","54","60")),], aes(x = 7, xend = 11, y = Angle.Fwd, yend = Angle.Fwd), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  geom_segment(data = hline_df[which(hline_df$QueryPoint %in% c("36","42","48","54","60")),], aes(x = 7, xend = 11, y = Angle.Back, yend = Angle.Back), 
               color = "blue", linetype = "dashed", size = 0.7, alpha=0.5) +
  geom_point(data = hline_df, aes(x = 1:11, y = Angle.Rev), color = "purple", size = 4, alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:11, y = Angle.Rev.180), color = "purple", size = 4, alpha=0.8) +
  theme_sticcc() +
  labs(x = "Trajectory Point", y = "Angle")
image

pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill.pdf")), width = 10, height = 10)
print(image)
dev.off()

# another plot to cut the legend off from
ggplot() +
  geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y,
                                   color = violinwidth), show.legend = FALSE) +
  # Re-use geom_violin to plot the outline
  geom_violin(data = rs_boxplot, aes(x = as.integer(QueryPoint), y = Angle, fill = QueryPoint),
              color = "white", alpha = 0, draw_quantiles = c(0.25, 0.5, 0.75),
              show.legend = TRUE) +
  scale_x_continuous(breaks = breaks, labels = labels) +
  scale_color_viridis_c(name="Density") 


# Density plot with reversibility
revScalingFactor <- 0.4
image <- ggplot() +
  #geom_point(data=pca$x, aes(x=PC1, y=PC2)) +
  geom_density2d(data=traj_pca[vector_subset,], aes(x=PC1, y=PC2), color="red") +
  geom_point(data=clust_431_traj, aes(x=PC1, y=PC2, color=1:60)) +
  geom_segment(data = traj_v_pred_REV, 
               aes(x=x,y=y, xend=x+dx*revScalingFactor, yend=y+dy*revScalingFactor), 
               arrow = arrow(length = unit(0.3,"cm")), color="black", size=2, alpha=0.7) +
  scale_color_gradient(name="Query Point", breaks=c(10,30,50)) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))


pdf(file = file.path(plotDir, paste0("VObs_On_PCA.pdf")), width = 10, height = 10)
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

#debug(v_obs_along_path)
rs_list <- v_obs_along_path(trajectory = traj_pca, # PCA plus Time column
                            lag = optimal_lags$OptimalLag, # numeric - time gap (computed numerically, so should correspond to times, not indices!)
                            sce = stic,
                            queryTrajectory = subset_trajectory, # either rownames of trajectory, or a dataframe of same ncol to be compared to it
                            neighborhoodRadius = 0.03,
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
selectedPtIdx <- 18
selectedPt <-  clust_431_traj[selectedPtIdx,]

image <- plot_v_vs_trajectory(sce = stic,
                              queryPoint = selectedPt, 
                              neighborhoodRadius = 0.02,
                              trajectory = traj_pca,
                              lag = 0.6,
                              v_pred = traj_v_pred_REV[4,], 
                              scale_factor = 1)
image




