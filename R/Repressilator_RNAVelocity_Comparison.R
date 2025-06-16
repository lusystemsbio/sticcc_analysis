##### GLOBALS & SETUP #####
rm(list=ls())
library(sRACIPE)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(STICCC)
library(dplyr)
library(tibble)
library(FNN)  # For k-nearest neighbors
set.seed(123)
source("R/sticcc_analysis_utilities.R")


# global params
topoName <- "repressilator"
forceSim <- FALSE     
forcePCA <- FALSE
forceSTICCC <- FALSE
saveNetworkPlot <- FALSE
nSamples <- 10000
pseudocount <- T
numClusters <- 6
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# directory setup
expName <- paste0(topoName,"_RNAVelocity")
topoDir <- file.path(getwd(),topoName)
expDir <- file.path(topoDir,"RNAVelocity_comparison")
pythonDataDir <- "/Users/danramirez/Desktop/NEU/projects/sticcc/python_export_data"

outputDir = file.path(expDir,"data")
plotDir <- file.path(expDir, "plots")

if(!dir.exists(topoDir)) {
  dir.create(topoDir)
}
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}
if(!dir.exists(plotDir)) {
  dir.create(plotDir)
}
if(!dir.exists(expDir)) {
  dir.create(expDir)
}

# load topology file
topo <- loadTopo(topoName)
topo$Type[which(topo$Type %% 2 == 0)] = 2
topo$Type[which(topo$Type %% 2 == 1)] = 1


##### IMPORT SIMULATION DATA #####
# With the repressilator simulations, cell 1 happened to have identical counts. It will be removed to avoid
# breaking the distance calculation
keep_idx <- c(2:1000)

exprMat <- t(read.csv(file.path(pythonDataDir, paste0(topoName,"_spliced_counts_raw.csv")), 
                      row.names = 1))#[,keep_idx]
exprMat_norm <- t(read.csv(file.path(pythonDataDir, paste0(topoName,"_spliced_counts_norm.csv")), 
                           row.names = 1))#[,keep_idx]
timepoints_subset <- read.csv(file.path(pythonDataDir, paste0(topoName,"_5kSubset_times.csv")), 
                           row.names = 1)

## Ported from scVelo
# Interesting note: when you specify normalize_layers (or other kwargs) in scv.pp.filter_and_normalize, 
# it overrides the default value log=True and skips log1p normalization
# this seems insane!
normalize_per_cell_r <- function(exprMat, counts_per_cell_after = NULL) {
  # Compute total counts per cell
  counts_per_cell <- colSums(exprMat)
  
  # Determine scaling factor
  if (is.null(counts_per_cell_after)) {
    counts_per_cell_after <- median(counts_per_cell)  # Use median total UMI count
  }
  
  # Avoid division by zero
  counts_per_cell <- counts_per_cell / counts_per_cell_after
  counts_per_cell[counts_per_cell == 0] <- 1
  
  # Normalize each cell
  exprMat_norm <- sweep(exprMat, 2, counts_per_cell, FUN = "/")
  
  return(exprMat_norm)
}


##### IMPORT RNA VELOCITY #####
velocities <- read.csv(file.path(pythonDataDir, paste0(topoName,"_velocity_vectors.csv")), 
                         row.names = 1)#[keep_idx,]

##### PCA & STICCC SETUP #####
# create SCE object
## TODO: make the lines below into a small wrapper method createVIC()
stic <- sticSE(topo = topo, exprMat = exprMat, normData = exprMat_norm,
               topoName = topoName, expName = expName)


# add metadata
#undebug(prepMetadata)
stic <- prepMetadata(stic, exprMat_norm, cluster = T, k = numClusters)

# run PCA
stic <- runPCA(stic, save=T, overwrite=forcePCA, fname=file.path(outputDir,"PCA_res.Rds"))
stic@metadata$pca_data$rotation[,2] <- -1*stic@metadata$pca_data$rotation[,2]
reducedDim(stic,"PCA")[,2] <- -1*reducedDim(stic,"PCA")[,2]
colData(stic)$Time <- timepoints_subset$Time

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

stic@metadata$params$plotDim <- "PCA"



##### CELL-WISE COMPARISON #####

# PCA & plotting params
pca_df <- reducedDim(stic, "PCA")
pca_info <- stic@metadata$pca_data

scalingFactor <- 0.1

xMin <- stic@metadata$params$xMin
xMax <- stic@metadata$params$xMax
yMin <- stic@metadata$params$yMin
yMax <- stic@metadata$params$yMax

pc1_weight <- round(100*stic@metadata$pca_summary$importance[2,1],2)
pc2_weight <- round(100*stic@metadata$pca_summary$importance[2,2],2)
plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")

## First, plot each with individual vectors on PCA
st_vectors_v1 <- stic@metadata$vectors
st_vectors_v2 <- stic@metadata$vectors_in
st_net <- as.data.frame(st_vectors_v1 + st_vectors_v2)
st_rev <- as.data.frame((st_vectors_v1 - st_vectors_v2) / 2)


# Plot STICCC vectors
net_df <- merge(pca_df, st_net, by="row.names")
image <- ggplot(net_df, aes(x=PC1, y=PC2)) +
  geom_point(color="grey", alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("STICCC Net Flow on ",expName," cells")) +
  guides(alpha="none", size="none", color="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))  +
  geom_segment(aes(xend=PC1+dPC1*scalingFactor, yend=PC2+dPC2*scalingFactor),
               arrow = arrow(length = unit(0.2,"cm")))
image  

pdf(file = file.path(plotDir, paste0("pca_netflow.pdf")), width = 10, height = 10)
print(image)
dev.off()


rev_df <- merge(pca_df, st_rev, by="row.names")
image <- ggplot(rev_df, aes(x=PC1, y=PC2)) +
  geom_point(color="grey", alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("STICCC Reversibility on ",expName," cells")) +
  guides(alpha="none", size="none", color="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))  +
  geom_segment(aes(xend=PC1+dPC1*scalingFactor, yend=PC2+dPC2*scalingFactor),
               arrow = arrow(length = unit(0.2,"cm")))
image

pdf(file = file.path(plotDir, paste0("pca_rev.pdf")), width = 10, height = 10)
print(image)
dev.off()


# Plot RNA velocity
scvelo_scalingFactor <- 0.1
velocities_pca <- as.data.frame(as.matrix(velocities) %*% pca_info$rotation)
colnames(velocities_pca) <- paste0("d",colnames(velocities_pca))
scvelo_df <- merge(pca_df, velocities_pca, by="row.names")
image <- ggplot(scvelo_df, aes(x=PC1, y=PC2)) +
  geom_point(color="grey", alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("RNA velocity on ",expName," cells")) +
  guides(alpha="none", size="none", color="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))  +
  geom_segment(aes(xend=PC1+dPC1*scvelo_scalingFactor, yend=PC2+dPC2*scvelo_scalingFactor),
               arrow = arrow(length = unit(0.2,"cm")))
image

pdf(file = file.path(plotDir, paste0("pca_scvelo.pdf")), width = 10, height = 10)
print(image)
dev.off()


# Function to compute cosine similarity
cosine_similarity <- function(v1, v2) {
  if (any(is.na(v1)) || any(is.na(v2))) return(NA) # Handle missing values
  sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
}

# Function to compute dot product
dot_product <- function(v1, v2) {
  if (any(is.na(v1)) || any(is.na(v2))) return(NA) # Handle missing values
  sum(v1 * v2)
}

# Merge all dataframes by rownames (keeping all rows)
velocities_pca <- velocities_pca %>% rownames_to_column("cell_id")
st_net <- st_net %>% rownames_to_column("cell_id")
st_rev <- st_rev %>% rownames_to_column("cell_id")
all_vectors <- Reduce(function(x, y) full_join(x, y, by = "cell_id"), 
                   list(velocities_pca, st_net, st_rev)) %>%
  column_to_rownames(var = "cell_id")

# Extract column indices
pca_cols <- c("dPC1.x", "dPC2.x", "dPC3.x")   # scVelo vectors
net_cols <- c("dPC1.y", "dPC2.y", "dPC3.y")   # st_net vectors
rev_cols <- c("dPC1", "dPC2", "dPC3")         # st_rev vectors (last join keeps original names)

# Compute cosine similarity and dot product
vector_diffs <- all_vectors %>%
  rowwise() %>%
  mutate(
    cos_sim_net = cosine_similarity(c_across(all_of(pca_cols)), c_across(all_of(net_cols))),
    dot_prod_net = dot_product(c_across(all_of(pca_cols)), c_across(all_of(net_cols))),
    cos_sim_rev = cosine_similarity(c_across(all_of(pca_cols)), c_across(all_of(rev_cols))),
    dot_prod_rev = dot_product(c_across(all_of(pca_cols)), c_across(all_of(rev_cols)))
  ) %>%
  select(cos_sim_net, dot_prod_net, cos_sim_rev, dot_prod_rev)  # Keep only the final results
vector_diffs <- as.data.frame(vector_diffs)
rownames(vector_diffs) <- rownames(all_vectors)


# Plot distribution of cosine similarity & dot product
# Convert data to long format for ggplot
cos_sim_long <- reshape2::melt(vector_diffs, measure.vars = c("cos_sim_net", "cos_sim_rev"),
                                    variable.name = "Method", value.name = "Cosine Similarity")

# Rename legend entries
cos_sim_long$Method <- factor(cos_sim_long$Method, 
                                   levels = c("cos_sim_net", "cos_sim_rev"), 
                                   labels = c("Net", "Rev"))

# Plot histogram
image <- ggplot(cos_sim_long, aes(x = `Cosine Similarity`, fill = Method)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  scale_fill_manual(values = c("blue", "red")) +  # Custom colors
  labs(title = "Cosine Similarity, RNA Velocity vs STICCC", fill = "Method") +
  ylab("Count") +
  theme_minimal() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20))
image

pdf(file = file.path(plotDir, paste0("comparison_hist_cosine_similarity.pdf")), width = 10, height = 10)
print(image)
dev.off()

# plot on PCA
net_df_cos <- merge(pca_df, vector_diffs, by="row.names")
image <- ggplot(net_df_cos, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=cos_sim_net), alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("Cosine Similarity on ",expName," cells")) +
  labs(color="Cosine Similarity") +
  guides(alpha="none", size="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20)) 
image

pdf(file = file.path(plotDir, paste0("comparison_pca_cosine_similarity_netflow.pdf")), width = 10, height = 10)
print(image)
dev.off()

image <- ggplot(net_df_cos, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=dot_prod_net), alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("Dot Product on ",expName," cells")) +
  labs(color="Dot Product") +
  guides(alpha="none", size="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20)) 
image

pdf(file = file.path(plotDir, paste0("comparison_pca_dotprod_netflow.pdf")), width = 10, height = 10)
print(image)
dev.off()




dot_prod_long <- reshape2::melt(vector_diffs, measure.vars = c("dot_prod_net", "dot_prod_rev"),
                               variable.name = "Method", value.name = "Dot Product")

# Rename legend entries
dot_prod_long$Method <- factor(dot_prod_long$Method, 
                              levels = c("dot_prod_net", "dot_prod_rev"), 
                              labels = c("Net", "Rev"))

# Plot histogram
image <- ggplot(dot_prod_long, aes(x = `Dot Product`, fill = Method)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  scale_fill_manual(values = c("blue", "red")) +  # Custom colors
  labs(title = "Dot Product, RNA Velocity vs STICCC", fill = "Method") +
  ylab("Count") +
  xlim(-2, 2) + # this omits some extreme outliers (likely the points from the center)
  theme_minimal() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20))
image

pdf(file = file.path(plotDir, paste0("comparison_hist_dotprod.pdf")), width = 10, height = 10)
print(image)
dev.off()





##### LOCAL INCONSISTENCY #####

compute_local_inconsistency <- function(velocity_vectors, dist_matrix, k = 30) {
  # Ensure velocity vectors exist
  if (is.null(velocity_vectors)) stop("Velocity vectors are missing.")
  
  # Compute k-nearest neighbors (excluding self)
  knn_indices <- apply(dist_matrix, 1, order)[2:(k + 1), ]  # Get indices of k nearest neighbors
  
  # Initialize local inconsistency vector
  local_inconsistency <- numeric(nrow(velocity_vectors))
  
  # Compute local inconsistency for each cell
  for (i in 1:nrow(velocity_vectors)) {
    v_i <- velocity_vectors[i, , drop = FALSE]  # Current cell velocity vector
    v_neighbors <- velocity_vectors[knn_indices[, i], , drop = FALSE]  # Neighbors' velocity vectors
    
    # Compute cosine distances
    similarities <- apply(v_neighbors, 1, function(v_n) cosine_similarity(v_i, v_n))
    distances <- 1 - similarities  # Convert similarity to distance
    
    # Average cosine distance
    local_inconsistency[i] <- mean(distances, na.rm = TRUE)
  }
  
  return(local_inconsistency)
}

dist_matrix <- as.matrix(stats::dist(pca_df[,c("PC1","PC2","PC3")], method = "euclidean"))


force_local_inconsistency <- T
velo_pca_fname <- file.path(outputDir,"scVelo_pca_vectors.Rds")
if(!file.exists(velo_pca_fname) | force_local_inconsistency) {
  scVelo_inconsistency <- compute_local_inconsistency(velocities_pca[,c("dPC1","dPC2","dPC3")], 
                                                      dist_matrix)
  velocities_pca$Inconsistency <- scVelo_inconsistency
  saveRDS(velocities_pca, velo_pca_fname)
} else {
  velocities_pca <- readRDS(velo_pca_fname)
}

st_net_fname <- file.path(outputDir,"st_vectors_net.Rds")
if(!file.exists(st_net_fname) | force_local_inconsistency) {
  sticcc_inconsistency_net <- compute_local_inconsistency(st_net[,c("dPC1","dPC2","dPC3")], 
                                                          dist_matrix[st_net$cell_id, st_net$cell_id])
  st_net$Inconsistency <- sticcc_inconsistency_net  
  saveRDS(st_net, st_net_fname)
} else {
  st_net <- readRDS(st_net_fname)
}

st_rev_fname <- file.path(outputDir,"st_vectors_rev.Rds")
if(!file.exists(st_rev_fname) | force_local_inconsistency) {
  sticcc_inconsistency_rev <- compute_local_inconsistency(st_rev[,c("dPC1","dPC2","dPC3")], 
                                                          dist_matrix[st_net$cell_id, st_net$cell_id])
  st_rev$Inconsistency <- sticcc_inconsistency_rev
  saveRDS(st_rev, st_rev_fname)
} else {
  st_rev <- readRDS(st_rev_fname)
}






summary(st_net$Inconsistency)
summary(st_rev$Inconsistency)
summary(velocities_pca$Inconsistency)



image <- ggplot(pca_df[which(rownames(pca_df) %in% velocities_pca$cell_id),], aes(x=PC1, y=PC2)) +
  geom_point(aes(color=velocities_pca$Inconsistency), alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("scVelo Local Inconsistency on ",expName," cells")) +
  labs(color="Local Inconsistency") +
  guides(alpha="none", size="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20)) 
image

pdf(file = file.path(plotDir, paste0("pca_local_inconsistency_scvelo.pdf")), width = 10, height = 10)
print(image)
dev.off()



image <- ggplot(pca_df[st_net$cell_id,], aes(x=PC1, y=PC2)) +
  geom_point(aes(color=st_net$Inconsistency), alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("Net Flow Local Inconsistency on ",expName," cells")) +
  labs(color="Local Inconsistency") +
  guides(alpha="none", size="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20)) 
image

pdf(file = file.path(plotDir, paste0("pca_local_inconsistency_netflow.pdf")), width = 10, height = 10)
print(image)
dev.off()


image <- ggplot(pca_df[st_rev$cell_id,], aes(x=PC1, y=PC2)) +
  geom_point(aes(color=st_rev$Inconsistency), alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("Reversibility Local Inconsistency on ",expName," cells")) +
  labs(color="Local Inconsistency") +
  guides(alpha="none", size="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20)) 
image

pdf(file = file.path(plotDir, paste0("pca_local_inconsistency_rev.pdf")), width = 10, height = 10)
print(image)
dev.off()

##### TRAJECTORY VALIDATION #####

###### IMPORT TRAJECTORY ######
traj_raw <- t(read.csv(file.path(pythonDataDir, paste0(topoName,"_traj_raw.csv")), 
                      row.names = 1))
traj_norm <- log1p(normalize_per_cell_r(traj_raw))
#traj_pca <- as.data.frame(t(traj_norm) %*% pca_info$rotation)
traj_pca <- as.data.frame(scale(t(traj_norm), pca_info$center, pca_info$scale) %*% pca_info$rotation)

traj_times <- read.csv(file.path(pythonDataDir, paste0(topoName,"_traj_times.csv")))
traj_pca$Time <- traj_times$Time #round(seq(0, 5000, length.out=nrow(traj_pca)), 2)


#test_pca <- as.data.frame(t(exprMat_norm) %*% pca_info$rotation)
#test_pca_2 <- as.data.frame(scale(t(exprMat_norm), pca_info$center, pca_info$scale) %*% pca_info$rotation)

traj_raw_det <- t(read.csv(file.path(pythonDataDir, paste0(topoName,"_traj_raw_deterministic.csv")), 
                           row.names = 1))
traj_raw_time <- traj_raw_det[1,]
traj_raw_det <- traj_raw_det[-1,]
traj_norm_det <- log1p(normalize_per_cell_r(traj_raw_det))
#traj_pca_det <- as.data.frame(t(traj_norm_det) %*% pca_info$rotation)
traj_pca_det <- as.data.frame(scale(t(traj_norm_det), pca_info$center, pca_info$scale) %*% pca_info$rotation)
traj_pca_det$Time <- traj_raw_time
  


###### DETERMINISTIC LIMIT CYCLE ######
# We will compare predictions over the course of one limit cycle
# First, we need to find the approximate length of a limit cycle & extract corresponding data
# We will omit the first and last 1/4 of the simulation, then find peaks
library(pracma)
# peak_A_indices <- findpeaks(traj_norm_det[1,(ncol(traj_norm_det)/4):ncol(traj_norm_det)], 
#           nups = 5, ndowns = 1, minpeakheight = median(traj_norm_det[1,]) + sd(traj_norm_det[1,]))[,2]
min_PC2_indices <- findpeaks(-1*traj_pca_det[(nrow(traj_pca_det)/4):nrow(traj_pca_det),2], 
                            nups = 5, ndowns = 1, minpeakheight = median(traj_pca_det[,2]) + sd(traj_pca_det[,2]))[,2]


start_idx <- min_PC2_indices[1] + round(nrow(traj_pca_det)/4)
end_idx <- min_PC2_indices[2] + round(nrow(traj_pca_det)/4) - 1 # end a little before the cycle completes
t_step <- 210
  
subset_trajectory <- traj_pca_det[seq(start_idx,end_idx,t_step),]
rownames(subset_trajectory) <- rownames(traj_pca_det)[seq(start_idx,end_idx,t_step)]
subset_traj_placeholder <- subset_trajectory

# sanity check plot
ggplot(traj_pca[1:10000,], aes(x=PC1,y=PC2)) +
  geom_point() +
  geom_point(data=subset_trajectory[1,], aes(x=PC1, y=PC2), size=5, color="blue") +
  geom_point(data=subset_trajectory, aes(x=PC1,y=PC2), color="red", size=3)

###### IDENTIFY COMPARISON PTS ######
## In order to find a suitable comparison to RNA velocity, we'll select points from the noisy data
# which was used to calculate RNA velocity which are closest to those sampled from the deterministic trajectory

subset_trajectory <- subset_traj_placeholder


subset_trajectory_sub <- data.frame()
new_tps_list <- list()
new_idx_list <- list()
for(pt in 1:nrow(subset_trajectory)) {
  query <- subset_trajectory[pt,which(!colnames(subset_trajectory) == "Time")]
  nn <- FNN::get.knnx(pca_df, query, k=1)
  
  newPt <- pca_df[nn$nn.index,]
  newPtTime <- colnames(traj_norm)[nn$nn.index]
  new_tps_list[[pt]] <- newPtTime
  new_idx_list[[pt]] <- rownames(pca_df)[nn$nn.index]
  
  subset_trajectory_sub <- rbind(subset_trajectory_sub, newPt)
  
}
colnames(subset_trajectory_sub) <- c("PC1", "PC2", "PC3")
subset_trajectory_sub$Index <- 1:nrow(subset_trajectory_sub)
subset_trajectory_sub$Time <- as.numeric(unlist(new_tps_list))
subset_trajectory_sub$CellName <- unlist(new_idx_list)


subset_trajectory <- subset_trajectory_sub




###### STICC ON TRAJECTORY ######
# Calculate vectors for each subsetted point
stic@metadata$params$plotDim <- "PCA"
stic@metadata$params$nComponents <- 3 

#undebug(trajectorySmoothVectors)
traj_v_pred <- trajectorySmoothVectors(trajectory = subset_trajectory[,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))], # trajectory in PCA coordinates
                                       sce = stic,
                                       neighborhoodRadius = 0.05,
                                       invertV2 = T,
                                       vec.use = "net") 


traj_v_pred_REV <- trajectorySmoothVectors(trajectory = subset_trajectory[,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))], # trajectory in PCA coordinates
                                           sce = stic,
                                           neighborhoodRadius = 0.05,
                                           invertV2 = T,
                                           vec.use = "rev") 


###### SET LAG TIME ######
## Identify lag for each point to keep variance constant
var_traj_df <- data.frame()
#optimize_lags <- seq(0.1, 3, 1)
optimize_lags <- c(seq(0.2, 2.8, 0.2), seq(3, 40, 3))

for(subsetPtIdx in rownames(subset_trajectory)) {
  subsetPt <- subset_trajectory[subsetPtIdx,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))]
  
  #debug(vobs_var_by_t)
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


target_rmsd <- max(var_traj_df[which(var_traj_df$Lag == optimize_lags[1]), "RMSD"])

optimal_lags <- data.frame(QueryPoint = unique(var_traj_df$QueryPoint), OptimalLag = NA, Var = NA)
for(pt in unique(var_traj_df$QueryPoint)) {
  # Find lag where var is closest to target_var
  pt_vars <- var_traj_df[which(var_traj_df$QueryPoint == pt & var_traj_df$Lag <= 10),]
  pt_lag <- pt_vars[which.min(abs(pt_vars$RMSD - target_rmsd)), "Lag"]
  pt_minVar <- pt_vars[which.min(abs(pt_vars$RMSD - target_rmsd)), "RMSD"]
  
  optimal_lags[which(optimal_lags$QueryPoint == pt), "OptimalLag"] <- pt_lag
  optimal_lags[which(optimal_lags$QueryPoint == pt), "RMSD"] <- pt_minVar
}



ggplot(var_traj_df[which(var_traj_df$Lag < 2),]) +
  geom_point(aes(x=Lag, y=RMSD, color=QueryPoint), size=3) +
  geom_line(aes(x=Lag, y=RMSD, color=QueryPoint)) +
  theme_sticcc() +
  geom_hline(aes(yintercept=target_rmsd), color="red", lty="dashed", size=2) +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))



###### CALCULATE V_OBS ######
undebug(v_obs_along_path)
rs_list <- v_obs_along_path(trajectory = traj_pca, # PCA plus Time column
                            lag = optimal_lags$OptimalLag, # numeric - time gap (computed numerically, so should correspond to times, not indices!)
                            sce = stic,
                            queryTrajectory = subset_trajectory[,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))], # either rownames of trajectory, or a dataframe of same ncol to be compared to it
                            neighborhoodRadius = 0.03,
                            v_pred = traj_v_pred) 



rs_summary <- rs_list$Summary
rs_boxplot <- rs_list$BoxplotData
rs_boxplot$QueryPoint <- factor(rs_boxplot$QueryPoint, 
                                levels=as.character(sort(unique(as.numeric(rs_boxplot$QueryPoint)))))


###### PLOT V_OBS ######
# Scale angles and generate vectors
vObs_scaling_factor <- 0.5
rs_summary$Obs.Mag.1.Scaled <- rs_summary$Obs.Mag.1 / max(rs_summary$Obs.Mag.1) * vObs_scaling_factor
rs_summary$Obs.Mag.2.Scaled <- rs_summary$Obs.Mag.2 / max(rs_summary$Obs.Mag.2, na.rm = T) * vObs_scaling_factor

rs_summary$Obs.Vector.X.1 <- cos(rs_summary$Obs.Angle.1) * rs_summary$Obs.Mag.1.Scaled
rs_summary$Obs.Vector.Y.1 <- sin(rs_summary$Obs.Angle.1) * rs_summary$Obs.Mag.1.Scaled

rs_summary$Obs.Vector.X.2 <- cos(rs_summary$Obs.Angle.2) * rs_summary$Obs.Mag.2.Scaled
rs_summary$Obs.Vector.Y.2 <- sin(rs_summary$Obs.Angle.2) * rs_summary$Obs.Mag.2.Scaled


# Get scVelo vectors for selected points
subset_traj_rnaVelocity <- velocities_pca[which(velocities_pca$cell_id %in% subset_trajectory$CellName),]
subset_traj_rnaVelocity <- merge(subset_traj_rnaVelocity, pca_df, by.x="cell_id", by.y="row.names")
subset_traj_rnaVelocity <- subset_traj_rnaVelocity[match(subset_trajectory$CellName, subset_traj_rnaVelocity$cell_id), ]



# Using sampled ideal path, calculate a fwd and back angle for each point (exclude first and last)
hline_df <- data.frame(QueryPoint=rs_summary$QueryPoint,
                       Angle.Fwd = NA,
                       Angle.Back = NA,
                       Angle.NetFlow = NA,
                       Angle.Rev = NA,
                       Angle.Rev.180 = NA,
                       Angle.scVelo = NA
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
  } else if(i == nrow(subset_trajectory)) {
    diff_fwd <- subset_trajectory[1,c("PC1","PC2")] - subset_trajectory[i,c("PC1","PC2")]
    angle_fwd <- angle_conversion(diff_fwd)
  }
  
  # Get corresponding RNA velocity
  
  
  # Store output
  hline_df[i,"Angle.Fwd"] <- angle_fwd
  hline_df[i,"Angle.Back"] <- angle_back
  hline_df[i,"Angle.NetFlow"] <- angle_conversion(traj_v_pred[i, c("dx","dy")]) 
  hline_df[i,"Angle.Rev"] <- angle_conversion(traj_v_pred_REV[i, c("dx","dy")]) 
  hline_df[i,"Angle.Rev.180"] <- angle_conversion(-1*traj_v_pred_REV[i, c("dx","dy")]) 
  hline_df[i,"Angle.scVelo"] <- angle_conversion(subset_traj_rnaVelocity[i,c("dPC1", "dPC2")]) 
  
  
  
}


hline_df_copy <- hline_df


blue_segment_df <- data.frame(x=1, xend=10, y=hline_df$Angle.Fwd[1], yend=hline_df$Angle.Fwd[10])


# Violin plot with annotations
image <- ggplot(data=rs_boxplot) +
  geom_violin(aes(x=QueryPoint, y=Angle)) +
  #geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.Fwd), color = "blue", size = 4, alpha=0.8) +
  geom_segment(data = blue_segment_df, aes(x = x, xend = xend, y = y, yend=yend), color = "blue", size = 0.7, linetype="dashed", alpha=0.8) +
  geom_point(data = hline_df, aes(x = QueryPoint, y = Angle.NetFlow), color = "purple", size = 4, alpha=0.8) +
  labs(x = "Trajectory Point") +
  theme_sticcc() +
  theme(axis.text.x = element_text(angle=90))
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

library(scales)
image <- ggplot() +
  geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y,
                                   color = violinwidth), show.legend = FALSE) +
  # Re-use geom_violin to plot the outline
  geom_violin(data = rs_boxplot, aes(x = as.integer(QueryPoint), y = Angle, fill = QueryPoint),
              color = "white", alpha = 0, draw_quantiles = c(0.25, 0.5, 0.75),
              show.legend = FALSE) +
  scale_x_continuous(breaks = breaks, labels = labels) +
  scale_fill_discrete(guide = "none") +
  scale_color_viridis_c(name="Density") +
  #geom_segment(data = blue_segment_df, aes(x = x, xend = xend, y = y, yend=yend), color = "blue", size = 0.7, linetype="dashed", alpha=0.8) +
  geom_path(data = hline_df, aes(x = as.numeric(QueryPoint), y = Angle.Fwd, 
                                 linetype = "Det. Limit Cycle"), 
            color = "blue", size = 0.7, alpha=0.8, show.legend = TRUE) +
  geom_point(data = hline_df, aes(x = 1:10, y = Angle.NetFlow, shape = "STICCC"), size = 4, color = "purple", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:10, y = Angle.scVelo, shape = "scVelo"), size = 4, color = "#56B4E9", alpha=0.8) +
  guides(
    shape = guide_legend(override.aes = list(color = c("#56B4E9", "purple"), size = 4, alpha = 0.8))
  ) +
  theme_sticcc() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_linetype_manual(name = "Curve", values = c("Det. Limit Cycle" = "dashed")) +
  scale_shape_manual(name = "Method", values = c("STICCC" = 16, "scVelo" = 17)) +  # Different point shapes
  labs(x = "Trajectory Point", y = "Angle")
image



pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill_comparison_v2.pdf")), width = 10, height = 10)
print(image)
dev.off()


pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill_comparison_v2_NOLEGEND.pdf")), 
    width = 8.5, height = 10)
print(image + guides(shape="none", color="none", linetype="none"))
dev.off()




# Density plot with vectors from both methods
densScalingFactor_sticcc <- 1
densScalingFactor_scVelo <- 0.4

# Add method labels for each dataframe
traj_v_pred$Method <- "STICCC"
subset_traj_rnaVelocity$Method <- "scVelo"

# Combine the data
arrows_df <- rbind(
  traj_v_pred %>% dplyr::select(x, y, dx, dy, Method),
  subset_traj_rnaVelocity %>% dplyr::select(PC1, PC2, dPC1, dPC2, Method) %>% 
    rename(x = PC1, y = PC2, dx = dPC1, dy = dPC2)
)

# Create the plot
image <- ggplot() +
  geom_density2d(data=traj_pca[which(traj_pca$Time %in% timepoints_subset$Time),], 
                 aes(x=PC1, y=PC2), color="red") +
  geom_point(data=traj_pca_det[seq(start_idx,end_idx,40),], 
             aes(x=PC1, y=PC2), color="black") +
  geom_segment(data = arrows_df, 
               aes(x=x, y=y, xend=x+dx*(Method=="STICCC")*densScalingFactor_sticcc + dx*(Method=="scVelo")*densScalingFactor_scVelo, 
                   yend=y+dy*(Method=="STICCC")*densScalingFactor_sticcc + dy*(Method=="scVelo")*densScalingFactor_scVelo, 
                   color=Method), 
               arrow = arrow(length = unit(0.3,"cm")), size=2, alpha=1) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_color_manual(name="Method", values=c("STICCC"="purple", "scVelo"="#56B4E9")) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))

image

pdf(file = file.path(plotDir, paste0("VPred_On_PCA_comparison_v2.pdf")), width = 10, height = 10)
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







##### DROPOUT ROBUSTNESS #####














