##### GLOBALS & SETUP #####
rm(list=ls())
library(sRACIPE)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(STICCC)
library(dplyr)
library(tidyr)
library(tibble)
library(FNN)  # For k-nearest neighbors
set.seed(123)
source("R/sticcc_analysis_utilities.R")


# global params
topoName <- "CTS"
forceSim <- FALSE     
forcePCA <- FALSE
forceSTICCC <- FALSE
saveNetworkPlot <- FALSE
nSamples <- 10000
pseudocount <- T
numClusters <- 3
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
if(!dir.exists(expDir)) {
  dir.create(expDir)
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

##### IMPORT SIMULATION DATA #####
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
velocities_dynamical <- read.csv(file.path(pythonDataDir, paste0(topoName,"_velocity_vectors_dyn.csv")), 
                                 row.names = 1)
velocities_veloVI <- read.csv(file.path(pythonDataDir, paste0(topoName,"_velocity_vectors_velovi.csv")), 
                              row.names = 1)


##### PCA & STICCC SETUP #####
# create SCE object
## TODO: make the lines below into a small wrapper method createVIC()
stic <- sticSE(topo = topo, exprMat = exprMat, normData = exprMat_norm,
               topoName = topoName, expName = expName)


# add metadata
stic <- prepMetadata(stic, exprMat_norm, cluster = T, k = numClusters)

# add PCA from previous simulation
pca <- readRDS(file.path(topoDir,"data",paste0("PCA_res.Rds")))
#stic <- runPCA(stic, save=T, overwrite=forcePCA, fname=file.path(outputDir,"PCA_res.Rds"))
pca_df <- predict(pca, as.data.frame(t(exprMat_norm)))
reducedDim(stic,"PCA") <- pca_df
stic@metadata$pca_data <- pca[1:4]
stic@metadata$pca_summary <- summary(pca)


# compute grid based on PCA
stic <- computeGrid(stic)

# compute pairwise distance between points
stic <- computeDist(stic)



##### RUN STICCC #####
# compute trajectories
stic_fname <- file.path(outputDir, paste0("stic_",topoName,"_2025.Rds"))
if(!file.exists(stic_fname) | forceSTICCC) {
  stic <- runSTICCC(stic, v2=T, invertV2=T)
  saveRDS(stic, stic_fname)
} else {
  stic <- readRDS(stic_fname)
}

colData(stic)$Time <- timepoints_subset$Time


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


# Plot RNA velocity
scvelo_dyn_scalingFactor <- 0.2
velocities_dyn_pca <- as.data.frame(as.matrix(velocities_dynamical) %*% pca_info$rotation)
colnames(velocities_dyn_pca) <- paste0("d",colnames(velocities_dyn_pca))
scvelo_dyn_df <- merge(pca_df, velocities_dyn_pca, by="row.names")
image <- ggplot(scvelo_dyn_df, aes(x=PC1, y=PC2)) +
  geom_point(color="grey", alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("RNA velocity on ",expName," cells (sc-dyn)")) +
  guides(alpha="none", size="none", color="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))  +
  geom_segment(aes(xend=PC1+dPC1*scvelo_dyn_scalingFactor, yend=PC2+dPC2*scvelo_dyn_scalingFactor),
               arrow = arrow(length = unit(0.2,"cm")))
image

pdf(file = file.path(plotDir, paste0("pca_scvelo_dyn.pdf")), width = 10, height = 10)
print(image)
dev.off()


# Plot RNA velocity
velovi_scalingFactor <- 0.5
velocities_velovi_pca <- as.data.frame(as.matrix(velocities_veloVI) %*% pca_info$rotation)
colnames(velocities_velovi_pca) <- paste0("d",colnames(velocities_velovi_pca))
velovi_df <- merge(pca_df, velocities_velovi_pca, by="row.names")
image <- ggplot(velovi_df, aes(x=PC1, y=PC2)) +
  geom_point(color="grey", alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("RNA velocity on ",expName," cells (veloVI)")) +
  guides(alpha="none", size="none", color="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))  +
  geom_segment(aes(xend=PC1+dPC1*velovi_scalingFactor, yend=PC2+dPC2*velovi_scalingFactor),
               arrow = arrow(length = unit(0.2,"cm")))
image

pdf(file = file.path(plotDir, paste0("pca_velovi.pdf")), width = 10, height = 10)
print(image)
dev.off()



### Plot grid points
stic <- computeGrid(stic, grid.length = 15)
stic <- computeGridVectors(stic, inVectors = F, combine = T, unitVectors = F, how="net")

plotGrid(sce = stic,
                  colorVar = NA,
                  plotLoadings = F,
                  plotSuffix = paste0("_sep25_net_grey"),
                  minMagnitude = 0.001,
                  scalingFactor = 2,
                  arrowheadSize = 0.5,
                  return = T
                  
)


stic <- computeGridVectors(stic, inVectors = F, combine = T, unitVectors = F, how="rev")

plotGrid(sce = stic,
                  colorVar = NA,
                  plotLoadings = F,
                  plotSuffix = paste0("_sep25_rev_grey"),
                  minMagnitude = 0.001,
                  scalingFactor = 2,
                  arrowheadSize = 0.5,
                  return = T
                  
)



## scvelo grid plot

stic_scv <- stic
colData(stic_scv)$dX = velocities_pca$dPC1
colData(stic_scv)$dY = velocities_pca$dPC2
stic_scv <- computeGridVectors(stic_scv, inVectors = F, combine = F, unitVectors = F, how=NA)

plotGrid(sce = stic_scv,
         colorVar = NA,
         plotLoadings = F,
         plotSuffix = paste0("_sep25_scvelo"),
         minMagnitude = 0.001,
         scalingFactor = 0.3,
         arrowheadSize = 0.5,
         return = T
         
)

## scvelo dynamical
colData(stic_scv)$dX = velocities_dyn_pca$dPC1
colData(stic_scv)$dY = velocities_dyn_pca$dPC2
stic_scv <- computeGridVectors(stic_scv, inVectors = F, combine = F, unitVectors = F, how=NA)

plotGrid(sce = stic_scv,
         colorVar = NA,
         plotLoadings = F,
         plotSuffix = paste0("_sep25_scvelo_dyn"),
         minMagnitude = 0.001,
         scalingFactor = 0.6,
         arrowheadSize = 0.5,
         return = T
         
)


# veloVI
colData(stic_scv)$dX = velocities_velovi_pca$dPC1
colData(stic_scv)$dY = velocities_velovi_pca$dPC2
stic_scv <- computeGridVectors(stic_scv, inVectors = F, combine = F, unitVectors = F, how=NA)

plotGrid(sce = stic_scv,
         colorVar = NA,
         plotLoadings = F,
         plotSuffix = paste0("_sep25_veloVI"),
         minMagnitude = 0.001,
         scalingFactor = 2,
         arrowheadSize = 0.5,
         return = T
         
)



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
velocities_dyn_pca <- velocities_dyn_pca %>% rownames_to_column("cell_id")
velocities_velovi_pca <- velocities_velovi_pca %>% rownames_to_column("cell_id")
st_net <- st_net %>% rownames_to_column("cell_id")
st_rev <- st_rev %>% rownames_to_column("cell_id")
all_vectors <- Reduce(function(x, y) full_join(x, y, by = "cell_id"), 
                      list(velocities_pca, st_net, st_rev, velocities_dyn_pca, velocities_velovi_pca)) %>%
  column_to_rownames(var = "cell_id")

# Extract column indices
pca_cols <- c("dPC1.x", "dPC2.x", "dPC3.x")   # scVelo vectors
net_cols <- c("dPC1.y", "dPC2.y", "dPC3.y")   # st_net vectors
rev_cols <- c("dPC1.x.x", "dPC2.x.x", "dPC3.x.x")         # st_rev vectors 
scvelo_dyn_cols <- c("dPC1.y.y", "dPC2.y.y", "dPC3.y.y")   # scVelo dynamical vectors
velovi_cols <- c("dPC1", "dPC2", "dPC3")   # veloVI vectors  (last join keeps original names)

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
  ggtitle(paste0("Cosine Similarity on ",expName," cells (Net)")) +
  labs(color="Cosine Similarity") +
  guides(alpha="none", size="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20)) 
image


pdf(file = file.path(plotDir, paste0("comparison_pca_cosine_similarity_netflow.pdf")), width = 10, height = 10)
print(image)
dev.off()

image <- ggplot(net_df_cos, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=cos_sim_rev), alpha=0.8) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("Cosine Similarity on ",expName," cells (Rev)")) +
  labs(color="Cosine Similarity") +
  guides(alpha="none", size="none") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20)) 
image

pdf(file = file.path(plotDir, paste0("comparison_pca_cosine_similarity_rev.pdf")), width = 10, height = 10)
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


##### COMPARE BETWEEN ALL METHODS #####

# plot ordering & setup
method_order <- c("STICCC Rev", "STICCC Net", "scVelo", "scVelo dynamical", "veloVI")
method_cols <- list(
  "STICCC Net"       = net_cols,
  "STICCC Rev"       = rev_cols,
  "scVelo"           = pca_cols,
  "scVelo dynamical" = scvelo_dyn_cols,
  "veloVI"           = velovi_cols
)

# Helper to compute average of a per-cell vector similarity across overlapping cells
.avg_pairwise <- function(colsA, colsB, sim_fun) {
  cols_needed <- c(colsA, colsB)
  cols_needed <- cols_needed[cols_needed %in% colnames(all_vectors)]
  if (length(cols_needed) < 6) return(NA_real_)  # missing vectors
  ok <- stats::complete.cases(all_vectors[, cols_needed, drop = FALSE])
  if (!any(ok)) return(NA_real_)
  A <- as.matrix(all_vectors[ok, colsA, drop = FALSE])
  B <- as.matrix(all_vectors[ok, colsB, drop = FALSE])
  vals <- vapply(seq_len(nrow(A)), function(i) sim_fun(A[i, ], B[i, ]), numeric(1))
  mean(vals, na.rm = TRUE)
}

# Build matrix
cosine_mat <- matrix(NA_real_, nrow = length(method_order), ncol = length(method_order),
                     dimnames = list(method_order, method_order))

for (i in seq_along(method_order)) {
  for (j in seq_along(method_order)) {
    mi <- method_order[i]
    mj <- method_order[j]
    cosine_mat[i, j] <- .avg_pairwise(method_cols[[mi]], method_cols[[mj]], cosine_similarity)
  }
}

# Save the raw matrices for later use
cosine_df <- as.data.frame(cosine_mat, check.names = FALSE)
write.csv(cosine_df, file = file.path(outputDir, "pairwise_cosine_similarity_matrix_allMethods.csv"), row.names = TRUE)

# Plot heatmaps (ggplot2)
cos_long <- reshape2::melt(cosine_mat, varnames = c("MethodA","MethodB"), value.name = "AvgCosine")
cos_long_upper <- cos_long %>%
  dplyr::mutate(
    i = match(MethodA, method_order),
    j = match(MethodB, method_order)
  ) %>%
  dplyr::filter(j >= i)

# Plot only upper triangle
p_cos_upper <- ggplot(cos_long_upper, aes(x = MethodB, y = MethodA, fill = AvgCosine)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", AvgCosine)), size = 4) +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Avg Cosine") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16))

# Display
print(p_cos_upper)

# Save PDFs
ggsave(file.path(plotDir, "pairwise_cosine_similarity_heatmap.pdf"), p_cos_upper, width = 8, height = 6)




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

dist_matrix <- as.matrix(stats::dist(pca_df[,c("PC1","PC2","PC3","PC4")], method = "euclidean"))


force_local_inconsistency <- F
velo_pca_fname <- file.path(outputDir,"scVelo_pca_vectors.Rds")
if(!file.exists(velo_pca_fname) | force_local_inconsistency) {
  scVelo_inconsistency <- compute_local_inconsistency(velocities_pca[,c("dPC1","dPC2","dPC3","dPC4")], 
                                                      dist_matrix)
  velocities_pca$Inconsistency <- scVelo_inconsistency
  saveRDS(velocities_pca, velo_pca_fname)
} else {
  velocities_pca <- readRDS(velo_pca_fname)
}

st_net_fname <- file.path(outputDir,"st_vectors_net.Rds")
if(!file.exists(st_net_fname) | force_local_inconsistency) {
  sticcc_inconsistency_net <- compute_local_inconsistency(st_net[,c("dPC1","dPC2","dPC3","dPC4")], 
                                                          dist_matrix[st_net$cell_id, st_net$cell_id])
  st_net$Inconsistency <- sticcc_inconsistency_net  
  saveRDS(st_net, st_net_fname)
} else {
  st_net <- readRDS(st_net_fname)
}

st_rev_fname <- file.path(outputDir,"st_vectors_rev.Rds")
if(!file.exists(st_rev_fname) | force_local_inconsistency) {
  sticcc_inconsistency_rev <- compute_local_inconsistency(st_rev[,c("dPC1","dPC2","dPC3","dPC4")], 
                                                          dist_matrix[st_net$cell_id, st_net$cell_id])
  st_rev$Inconsistency <- sticcc_inconsistency_rev
  saveRDS(st_rev, st_rev_fname)
} else {
  st_rev <- readRDS(st_rev_fname)
}




summary(st_net$Inconsistency)
summary(st_rev$Inconsistency)
summary(velocities_pca$Inconsistency)



ggplot(pca_df[which(rownames(pca_df) %in% velocities_pca$cell_id),], aes(x=PC1, y=PC2)) +
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



ggplot(pca_df[st_net$cell_id,], aes(x=PC1, y=PC2)) +
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

ggplot(pca_df[st_rev$cell_id,], aes(x=PC1, y=PC2)) +
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

##### TRAJECTORY VALIDATION #####

###### IMPORT TRAJECTORY ######
traj_raw <- t(read.csv(file.path(pythonDataDir, paste0(topoName,"_traj_raw.csv")), 
                       row.names = 1))
traj_norm <- log1p(normalize_per_cell_r(traj_raw))
traj_pca <- as.data.frame(scale(t(traj_norm), pca_info$center, pca_info$scale) %*% pca_info$rotation)

traj_times <- read.csv(file.path(pythonDataDir, paste0(topoName,"_traj_times.csv")))
traj_pca$Time <- traj_times$Time 




###### DEFINE LINEAR PATH ######
# Based on this plot, we'll define a path from 1 --> 2 --> 3
clusters <- as.factor(colData(stic)$Cluster)
ggplot(pca_df, aes(x=PC1,y=PC2)) +
  geom_point(aes(color=clusters))


clust1_centroid <- colMeans(pca_df[which(clusters == 1),])
clust2_centroid <- colMeans(pca_df[which(clusters == 2),])
clust3_centroid <- colMeans(pca_df[which(clusters == 3),])


ggplot(pca_df, aes(x=PC1,y=PC2)) +
  geom_point(aes(color=as.factor(clusters))) +
  geom_point(data=t(as.data.frame(clust1_centroid)), aes(x=PC1, y=PC2), size=3, color="black") +
  geom_point(data=t(as.data.frame(clust2_centroid)), aes(x=PC1, y=PC2), size=3, color="black") +
  geom_point(data=t(as.data.frame(clust3_centroid)), aes(x=PC1, y=PC2), size=3, color="black") 
  
  

clust_transit_traj <- distinct(rbind(interpolate_vectors(clust1_centroid, clust2_centroid, 5), 
                            interpolate_vectors(clust2_centroid, clust3_centroid, 5)))
colnames(clust_transit_traj) <- colnames(pca$x)


ggplot() +
  #geom_point(data=pca$x, aes(x=PC1, y=PC2)) +
  geom_density2d(data=pca_df, aes(x=PC1, y=PC2)) +
  geom_point(data=clust_transit_traj, aes(x=PC1, y=PC2, color=1:nrow(clust_transit_traj))) +
  scale_color_gradient(name="Query Point") +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))



###### IDENTIFY COMPARISON PTS ######
## In order to find a suitable comparison to RNA velocity, we'll select points from the noisy data
# which was used to calculate RNA velocity which are closest to those sampled from the deterministic trajectory

subset_trajectory <- clust_transit_traj


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
colnames(subset_trajectory_sub) <- c("PC1", "PC2", "PC3", "PC4")
subset_trajectory_sub$Index <- 1:nrow(subset_trajectory_sub)
subset_trajectory_sub$Time <- as.numeric(unlist(new_tps_list))
subset_trajectory_sub$CellName <- unlist(new_idx_list)


subset_trajectory <- subset_trajectory_sub

###### STICC ON TRAJECTORY ######
stic@metadata$params$plotDim <- "PCA"
stic@metadata$params$nComponents <- 4 
# Calculate vectors for each subsetted point
traj_v_pred <- trajectorySmoothVectors(trajectory = subset_trajectory[,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))], # trajectory in PCA coordinates
                                       sce = stic,
                                       neighborhoodRadius = 0.03,
                                       invertV2 = T,
                                       vec.use = "net") 


traj_v_pred_REV <- trajectorySmoothVectors(trajectory = subset_trajectory[,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))], # trajectory in PCA coordinates
                                           sce = stic,
                                           neighborhoodRadius = 0.03,
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
  geom_hline(aes(yintercept=target_rmsd), color="red", lty="dashed", linewidth=2) +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))


###### CALCULATE V_OBS ######
undebug(v_obs_along_path)
rs_list <- v_obs_along_path(trajectory = traj_pca, # PCA plus Time column
                            lag = optimal_lags$OptimalLag, # numeric - time gap (computed numerically, so should correspond to times, not indices!)
                            sce = stic,
                            queryTrajectory = subset_trajectory[,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))], # either rownames of trajectory, or a dataframe of same ncol to be compared to it
                            neighborhoodRadius = 0.01, # given the large trajectory size, we want a very small radius
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

subset_traj_scv_dyn <- velocities_dyn_pca[which(velocities_dyn_pca$cell_id %in% subset_trajectory$CellName),]
subset_traj_scv_dyn <- merge(subset_traj_scv_dyn, pca_df, by.x="cell_id", by.y="row.names")
subset_traj_scv_dyn <- subset_traj_scv_dyn[match(subset_trajectory$CellName, subset_traj_scv_dyn$cell_id), ]

subset_traj_velovi <- velocities_velovi_pca[which(velocities_velovi_pca$cell_id %in% subset_trajectory$CellName),]
subset_traj_velovi <- merge(subset_traj_velovi, pca_df, by.x="cell_id", by.y="row.names")
subset_traj_velovi <- subset_traj_velovi[match(subset_trajectory$CellName, subset_traj_velovi$cell_id), ]


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
    # this correction only applies for cyclic case, i.e. repressilator
    #diff_fwd <- subset_trajectory[1,c("PC1","PC2")] - subset_trajectory[i,c("PC1","PC2")]
    #angle_fwd <- angle_conversion(diff_fwd)
  }
  
  # Get corresponding RNA velocity
  # Store output
  hline_df[i,"Angle.Fwd"] <- angle_fwd
  hline_df[i,"Angle.Back"] <- angle_back
  hline_df[i,"Angle.NetFlow"] <- angle_conversion(traj_v_pred[i, c("dx","dy")]) 
  hline_df[i,"Angle.Rev"] <- angle_conversion(traj_v_pred_REV[i, c("dx","dy")]) 
  hline_df[i,"Angle.Rev.180"] <- angle_conversion(-1*traj_v_pred_REV[i, c("dx","dy")]) 
  hline_df[i,"Angle.scVelo"] <- angle_conversion(subset_traj_rnaVelocity[i,c("dPC1", "dPC2")]) 
  hline_df[i,"Angle.scVelo.Dyn"] <- angle_conversion(subset_traj_scv_dyn[i,c("dPC1", "dPC2")]) 
  hline_df[i,"Angle.veloVI"] <- angle_conversion(subset_traj_velovi[i,c("dPC1", "dPC2")]) 
  
  
  
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

# pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD.pdf")), width = 10, height = 10)
# print(image)
# dev.off()


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
              show.legend = TRUE) +
  scale_x_continuous(breaks = breaks, labels = labels) +
  scale_fill_discrete(guide = "none") +
  scale_color_viridis_c(name="Density") +
  # geom_path(data = hline_df, aes(x = as.numeric(QueryPoint), y = Angle.Fwd,
  #                                linetype = "Det. Limit Cycle"),
  #           color = "blue", size = 0.7, alpha=0.8, show.legend = TRUE, linetype="dashed") +
  #geom_path(data = hline_df, aes(x = as.numeric(QueryPoint), y = Angle.Back,
  #                               linetype = "Det. Limit Cycle"),
  #          color = "blue", size = 0.7, alpha=0.8, show.legend = TRUE, linetype="dashed") +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.Rev, shape = "STICCC"), size = 4, color = "purple", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.Rev.180, shape = "STICCC"), size = 4, color = "purple", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.scVelo, shape = "scVelo"), size = 4, color = "#56B4E9", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.scVelo.Dyn, shape = "scVelo-dyn"), size = 4, color = "#F0E442", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.veloVI, shape = "veloVI"), size = 4, color = "#D55E00", alpha=0.8) +
  guides(
    shape = guide_legend(override.aes = list(color = c("purple", "#56B4E9", "#F0E442", "#D55E00"), size = 4, alpha = 0.8))
  ) +
  theme_sticcc() +
  theme(axis.text.x = element_text(angle=90)) +
  #scale_linetype_manual(name = "Curve", values = c("Det. Limit Cycle" = "dashed")) +
  scale_shape_manual(name = "Method", values = c("STICCC" = 16, "scVelo" = 17, "scVelo-dyn"=4, "veloVI"=12)) +  # Different point shapes
  labs(x = "Trajectory Point", y = "Angle")
image



# pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill_comparison_v2_allMethods.pdf")), width = 10, height = 10)
# print(image)
# dev.off()
# 
# pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill_comparison_v2_NOLEGEND.pdf")), 
#     width = 8.5, height = 10)
# print(image + guides(shape="none", color="none", linetype="none"))
# dev.off()






# Density plot with vectors from both methods
densScalingFactor_sticcc <- 5
densScalingFactor_scVelo <- 1
densScalingFactor_scVelo_dyn <- 1
densScalingFactor_velovi <- 1

# Add method labels for each dataframe
traj_v_pred$Method <- "STICCC"
traj_v_pred_REV$Method <- "STICCC"
subset_traj_rnaVelocity$Method <- "scVelo"
subset_traj_scv_dyn$Method <- "scVelo-dyn"
subset_traj_velovi$Method <- "veloVI"

# Combine the data
arrows_df <- rbind(
  traj_v_pred_REV %>% dplyr::select(x, y, dx, dy, Method),
  subset_traj_rnaVelocity %>% dplyr::select(PC1, PC2, dPC1, dPC2, Method) %>% 
    rename(x = PC1, y = PC2, dx = dPC1, dy = dPC2),
  subset_traj_scv_dyn %>% dplyr::select(PC1, PC2, dPC1, dPC2, Method) %>% 
    rename(x = PC1, y = PC2, dx = dPC1, dy = dPC2),
  subset_traj_velovi %>% dplyr::select(PC1, PC2, dPC1, dPC2, Method) %>% 
    rename(x = PC1, y = PC2, dx = dPC1, dy = dPC2)
  
)
arrows_df$Method <- factor(arrows_df$Method, levels=c("STICCC","scVelo","scVelo-dyn","veloVI"))

image <- ggplot() +
  geom_density2d(data=traj_pca[which(traj_pca$Time %in% timepoints_subset$Time),], 
                 aes(x=PC1, y=PC2), color="red") +
  # STICCC arrows with ends="both"
  geom_segment(data = arrows_df[arrows_df$Method == "STICCC", ], 
               aes(x=x, y=y, 
                   xend=x+dx*densScalingFactor_sticcc, 
                   yend=y+dy*densScalingFactor_sticcc, 
                   color=Method), 
               arrow = arrow(length = unit(0.3,"cm"), ends = "both"), 
               size=2, alpha=0.8) +
  # scVelo arrows with ends="last" (or "first" if desired)
  geom_segment(data = arrows_df[arrows_df$Method != "STICCC", ], 
               aes(x=x, y=y, 
                   xend=x+dx*(Method=="scVelo")*densScalingFactor_scVelo +
                     dx*(Method=="scVelo-dyn")*densScalingFactor_scVelo_dyn + 
                     dx*(Method=="veloVI")*densScalingFactor_velovi, 
                   yend=y+dy*(Method=="scVelo")*densScalingFactor_scVelo +
                     dy*(Method=="scVelo-dyn")*densScalingFactor_scVelo_dyn +
                     dy*(Method=="veloVI")*densScalingFactor_velovi, 
                   color=Method), 
               arrow = arrow(length = unit(0.3,"cm"), ends = "last"),  # Change to "first" if needed
               size=2, alpha=0.8) +
  geom_point(data=subset_trajectory, 
             aes(x=PC1, y=PC2), color="black") +
  scale_color_manual(name="Method", values=c("STICCC"="purple", "scVelo"="#56B4E9", "scVelo-dyn"="#F0E442", "veloVI"="#D55E00")) +
  theme_sticcc() +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  theme(axis.line = element_line(linewidth = 0.7, colour = "black"))

image

pdf(file = file.path(plotDir, paste0("VPred_On_PCA_comparison_v2_allMethods.pdf")), width = 10, height = 10)
print(image)
dev.off()




ggplot(traj_pca[1:10000,]) +
  geom_point(data=pca_df, aes(x=PC1, y=PC2), color="gray") +
  geom_path(aes(x=PC1, y=PC2))




###### ACCURACY ACROSS METHODS ######

angdiff <- function(a, b) {
  # returns minimal absolute difference in radians in [0, pi]
  abs(atan2(sin(a - b), cos(a - b)))
}



## Find peak in v_obs angle distribution
## Calculate cosine similarity b/w peak angle and predicted vector
# If multiple predictions or multiple peaks exist, save highest score
# Compare against "theoretical range" of all possible predicted angles

library(circular)

#cos_sim  <- function(delta) cos(delta)
sim_agreement_smooth <- function(delta) (1 + cos(delta)) / 2


wrap_pi <- function(a) ((a + pi) %% (2*pi)) - pi




peak_results_df <- data.frame(QueryPoint=1:nrow(subset_trajectory),
                              Peak.Angle.1=NA,
                              Peak.Angle.2=NA,
                              Peak.Angle.3=NA,
                              n.Peaks=NA,
                              path_score=NA,
                              path_back_score=NA,
                              STICCC_net_score=NA,
                              STICCC_rev_score=NA,
                              scVelo_score=NA,
                              scVelo_dyn_score=NA,
                              veloVI_score=NA
)

theoretical_scores <- list()

for(subsetPtIdx in rownames(subset_trajectory)) {
  subsetPt <- subset_trajectory[subsetPtIdx,which(!colnames(subset_trajectory) %in% c("Time", "Index", "CellName"))]
  
  ## 1. Find peak angle of v_obs
  v_obs_dist <- rs_boxplot[which(rs_boxplot$QueryPoint == as.numeric(subsetPtIdx)),"Angle"]
  
  # circular KDE
  th <- circular(v_obs_dist, units = "radians", modulo = "2pi")
  dk <- density.circular(th, kernel = "vonmises", bw=4, n=2048)

  # Detect peaks
  peaks <- findpeaks(dk$y)

  # Get corresponding peak angles
  peak_angles <- dk$x[peaks[, 2]]
  peak_angles <- as.numeric( (peak_angles + pi) %% (2*pi) - pi )
  
  # put x on [-pi, pi] and keep x,y aligned & sorted
  x <- wrap_pi(as.numeric(dk$x))
  y <- as.numeric(dk$y)
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  
  # peak_angles <- peaks_from_diffs(x = x,
  #                                 y = y,
  #                                 prom_frac = 0.0000001,  
  #                                 min_sep   = 0.05) 

  d <- density(c(v_obs_dist),
               n = 4096, from = -pi, to = pi)
  plot(d$x, d$y, type = "l", xlab = "Angle", ylab = "Density")
  abline(v = peak_angles, lty = 2, col = "red")
  title(paste0("Point ",subsetPtIdx))
  
  
  plot(x, y, type = "l", xlab = "Angle", ylab = "Density")
  abline(v = peak_angles, lty = 2, col = "red")
  title(paste0("Point ",subsetPtIdx))
  
  
  
  # 2. Calculate best score
  # get predicted angle for 1) STICCC (net), 2) STICCC (rev), 3) STICCC (rev+180º), 4) scVelo, 5) scVelo-dyn, 6) veloVI
  pred_angles <- hline_df[subsetPtIdx,c(2:ncol(hline_df))]
  
  
  method_scores <- sapply(pred_angles, function(theta_pred) {
    theta_pred <- as.numeric(theta_pred)
    if (length(peak_angles) == 1L) {
      mean(sim_agreement_smooth(angdiff(theta_pred, peak_angles)), na.rm = TRUE)
    } else {
      # For each prediction, choose the best-matching peak
      diffs <- outer(theta_pred, peak_angles, function(t, p) angdiff(t, p))
      sims  <- sim_agreement_smooth(diffs)                         # rows: predictions, cols: peaks
      best_per_pred <- apply(sims, 1L, max, na.rm = TRUE)
      mean(best_per_pred, na.rm = TRUE)
    }
  })
  
  
  
  # 3. Calculate score range
  possible_angles <- seq(-2*pi, 2*pi, 0.02)
  
  possible_scores <- sapply(possible_angles, function(theta_pred) {
    theta_pred <- as.numeric(theta_pred)
    if (length(peak_angles) == 1L) {
      mean(sim_agreement_smooth(angdiff(theta_pred, peak_angles)), na.rm = TRUE)
    } else {
      # For each prediction, choose the best-matching peak
      diffs <- outer(theta_pred, peak_angles, function(t, p) angdiff(t, p))
      sims  <- sim_agreement_smooth(diffs)                         # rows: predictions, cols: peaks
      best_per_pred <- apply(sims, 1L, max, na.rm = TRUE)
      mean(best_per_pred, na.rm = TRUE)
    }
  })
  
  
  
  peak_results_df[subsetPtIdx,"n.Peaks"] <- length(peak_angles)
  peak_results_df[subsetPtIdx,"Peak.Angle.1"] <- peak_angles[1]
  if(length(peak_angles) > 1) {
    peak_results_df[subsetPtIdx,"Peak.Angle.2"] <- peak_angles[2]
    if(length(peak_angles) > 2) {
      peak_results_df[subsetPtIdx,"Peak.Angle.3"] <- peak_angles[3]
    }
  }
  
  peak_results_df[subsetPtIdx,"path_score"] <- method_scores["Angle.Fwd"]
  peak_results_df[subsetPtIdx,"path_back_score"] <- method_scores["Angle.Back"]
  peak_results_df[subsetPtIdx,"STICCC_net_score"] <- method_scores["Angle.NetFlow"]
  peak_results_df[subsetPtIdx,"STICCC_rev_score"] <- max(method_scores["Angle.Rev"],  method_scores["Angle.Rev.180"])
  peak_results_df[subsetPtIdx,"scVelo_score"] <- method_scores["Angle.scVelo"]
  peak_results_df[subsetPtIdx,"scVelo_dyn_score"] <- method_scores["Angle.scVelo.Dyn"]
  peak_results_df[subsetPtIdx,"veloVI_score"] <- method_scores["Angle.veloVI"]
  theoretical_scores[[subsetPtIdx]] <- possible_scores
  
  
  
}



# plot results
qp_names <- names(theoretical_scores)

theo_df <- do.call(
  rbind,
  Map(function(scores, qp)
    data.frame(QueryPoint = as.integer(qp), Score = as.numeric(scores)),
    theoretical_scores, qp_names)
)

theo_df_summary <- data.frame(QueryPoint=1:nrow(subset_trajectory),Max=NA,Min=NA)
for(i in 1:nrow(theo_df_summary)) {
  theo_df_summary[i,"Max"] <- max(theo_df[which(theo_df$QueryPoint == i),"Score"])
  theo_df_summary[i,"Min"] <- min(theo_df[which(theo_df$QueryPoint == i),"Score"])
}
theo_df_summary$QueryPoint <- factor(theo_df_summary$QueryPoint, levels = sort(unique(theo_df_summary$QueryPoint)))

# 2) Method scores long-form
method_cols <- c("STICCC_rev_score",
                 "scVelo_score","scVelo_dyn_score","veloVI_score")

scores_long <- peak_results_df %>%
  select(QueryPoint, all_of(method_cols)) %>%
  pivot_longer(-QueryPoint, names_to = "Method", values_to = "Score") %>%
  mutate(Method = factor(Method, levels = method_cols))

scores_long$Method <- factor(scores_long$Method, levels = c("STICCC_rev_score","scVelo_score","scVelo_dyn_score","veloVI_score"),
                             labels=c("STICCC (rev)", "scVelo", "scVelo-dyn", "veloVI"))

## ---------- Plot 1: per-QueryPoint boxplot + method points ----------
image <- ggplot() +
  geom_point(
    data = scores_long,
    aes(x = factor(QueryPoint),
        y = Score, color = Method, shape = Method, fill=Method),
    position = position_jitter(width = 0.15, height = 0, seed = 1),
    size = 6, alpha = 0.7
  ) +
  geom_rect(
    data = theo_df_summary,
    aes(xmin = as.integer(QueryPoint) - 0.45,
        xmax = as.integer(QueryPoint) + 0.45,
        ymin = Min, ymax = Max),
    color = "gray", alpha = 0.3
  ) +
  geom_point(
    data = scores_long,
    aes(x = factor(QueryPoint),
        y = Score, color = Method, shape = Method, fill=Method),
    position = position_jitter(width = 0.15, height = 0, seed = 1),
    size = 6, alpha = 0.7
  ) +
  labs(x = "Trajectory Point", y = "Score") +
  theme_sticcc() +
  scale_color_manual(values=c("purple", "#56B4E9", "#F0E442", "#D55E00")) +
  scale_fill_manual(values=c("purple", "#56B4E9", "#F0E442", "#D55E00")) +
  scale_shape_manual(name = "Method", values = c("STICCC (rev)" = 16, "scVelo" = 17, "scVelo-dyn"=25, "veloVI"=15), 
                     breaks=c("STICCC (rev)", "scVelo", "scVelo-dyn","veloVI")) +  # Different point shapes
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(shape="none")
image

pdf(file = file.path(plotDir, paste0("angle_scores_per_trajectoryPoint.pdf")), width = 10, height = 10)
print(image)
dev.off()



## ---------- Plot 2: per-Method boxplot over all QueryPoints ----------
library(ggpubr)

scores_long_boxplot <- scores_long %>% mutate(QueryPoint = as.factor(QueryPoint))


image <- ggplot(scores_long_boxplot, aes(x = Method, y = Score, fill = Method)) +
  geom_boxplot(outlier.shape = 16, width = 0.7) +
  scale_fill_manual(values=c("purple", "#56B4E9", "#F0E442", "#D55E00")) +
  theme_sticcc() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
image

pdf(file = file.path(plotDir, paste0("angle_scores_per_method.pdf")), width = 10, height = 10)
print(image)
dev.off()


pdf(file = file.path(plotDir, paste0("angle_scores_per_method_significance.pdf")), width = 10, height = 10)
print(image +
        geom_pwc(method = "wilcox_test",
                 p.adjust.method="none",
                 method.args = list(paired=TRUE),
                 label = "p.signif",
                 hide.ns = TRUE) +
        labs(x = "Method", y = "Score"))
dev.off()


## ----------- Plot 3: obs vs pred angle distributions -----------

mywidth <- 0.4
vl_fill_circ <- list()

query_levels <- sort(unique(rs_boxplot$QueryPoint))

for (i in seq_along(query_levels)) {
  qp <- query_levels[i]
  dsub <- subset(rs_boxplot, QueryPoint == qp)
  
  # circular KDE on [-pi, pi]
  th <- circular(dsub$Angle, units = "radians", modulo = "2pi")
  dens <- density.circular(th, kernel = "vonmises", bw=4, n = 512)
  
  
  y  <- as.numeric(dens$x)
  y[y > pi] <- y[y > pi] - 2*pi       # wrap if necessary
  d  <- as.numeric(dens$y)
  d  <- d / max(d)                    # normalize like violinwidth
  
  x     <- as.integer(i)
  xnew  <- x - mywidth * d
  xend  <- x + mywidth * d
  
  vl_fill_circ[[i]] <- tibble(
    QueryPoint = qp,
    x = x,
    y = y,
    violinwidth = d,
    xnew = xnew,
    xend = xend
  )
}

vl_fill_circ <- do.call(rbind, vl_fill_circ)


# Left and right “sides” of the violin, then stitch into a closed polygon per group
outline_list <- list()
ql <- sort(unique(vl_fill_circ$QueryPoint))

for (i in seq_along(ql)) {
  qp <- ql[i]
  d  <- vl_fill_circ[vl_fill_circ$QueryPoint == qp, c("y","xnew","xend")]
  o  <- order(d$y)                 # ensure monotonic y (angle)
  
  # left and right edges, then stitch into a closed polygon
  x_poly <- c(d$xnew[o], rev(d$xend[o]))
  y_poly <- c(d$y[o],    rev(d$y[o]))
  
  outline_list[[i]] <- data.frame(
    QueryPoint = qp,
    x_poly = x_poly,
    y_poly = y_poly
  )
}

outline_df <- do.call(rbind, outline_list)



quant_df <- rs_boxplot |>
  dplyr::group_by(QueryPoint) |>
  dplyr::reframe({
    th <- circular(Angle, type = "angles", units = "radians",
                   template = "none", modulo = "2pi")
    qs <- as.numeric(quantile.circular(th, probs = c(0.25, 0.5, 0.75)))
    qs[qs > pi] <- qs[qs > pi] - 2*pi       # wrap if necessary
    tibble::tibble(
      Quant = c("Q1","Q2","Q3"),
      Angle = qs,
      x     = as.integer(dplyr::first(QueryPoint))  # scalar!
    )
  }) |>
  dplyr::mutate(
    xnew = x - mywidth,
    xend = x + mywidth
  )



peak_results_df_long <- pivot_longer(peak_results_df[,c("QueryPoint","Peak.Angle.1","Peak.Angle.2","Peak.Angle.3")], 
                                     cols = c("Peak.Angle.1","Peak.Angle.2","Peak.Angle.3"), 
                                     values_to = "Angle",
                                     names_to = "PeakNum")







image <- ggplot() +
  geom_segment(data = vl_fill_circ, aes(x = xnew, xend = xend, y = y, yend = y,
                                        color = violinwidth), show.legend = FALSE) +
  # Re-use geom_violin to plot the outline
  geom_polygon(data = outline_df,
               aes(x = x_poly, y = y_poly, group = QueryPoint, fill = factor(QueryPoint)),
               color = "white", fill = NA, linewidth = 0.6, show.legend = FALSE) +
  scale_x_continuous(breaks = breaks, labels = labels) +
  geom_segment(data = quant_df,
               aes(x = xnew, xend = xend, y = Angle, yend = Angle),
               color = "white", linewidth = 0.6, alpha = 0.9, inherit.aes = FALSE) +
  scale_fill_discrete(guide = "none") +
  scale_color_viridis_c(name="Density") +
  geom_segment(
    data = peak_results_df_long,
    aes(x = as.integer(QueryPoint) - mywidth,
        xend = as.integer(QueryPoint) + mywidth,
        y = Angle, yend = Angle, linetype="Peak Angle"),
    color = "red", linewidth = 1.2, alpha = 0.8
  ) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.Rev, shape = "STICCC"), size = 4, color = "purple", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.Rev.180, shape = "STICCC"), size = 4, color = "purple", alpha=0.8) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.scVelo, shape = "scVelo"), size = 6, color = "#56B4E9", alpha=0.7) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.scVelo.Dyn, shape = "scVelo-dyn"), size = 6, color = "#F0E442", fill="#F0E442", alpha=0.7) +
  geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.veloVI, shape = "veloVI"), size = 6, color = "#D55E00", alpha=0.7) +
  guides(
    shape = guide_legend(override.aes = list(color = c("purple", "#56B4E9", "#F0E442",  "#D55E00"), size = 4, alpha = 0.8, breaks=c("STICCC", "scVelo", "scVelo-dyn","veloVI")))
  ) +
  theme_sticcc() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_linetype_manual(name = "Curve", values = c("Peak Angle" = "dotted")) +
  scale_shape_manual(name = "Method", values = c("STICCC" = 16, "scVelo" = 17, "scVelo-dyn"=25, "veloVI"=15), breaks=c("STICCC", "scVelo", "scVelo-dyn","veloVI")) +  # Different point shapes
  labs(x = "Trajectory Point", y = "Angle")
image


pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill_comparison_v2_allMethods_peaks.pdf")), width = 10, height = 10)
print(image)
dev.off()


pdf(file = file.path(plotDir, paste0("Angle_Dists_Optimized_Lag_RMSD_gradientFill_comparison_v2_NOLEGEND_allMethods_peaks.pdf")), 
    width = 10, height = 10)
print(image + guides(shape="none", color="none", linetype="none"))
dev.off()












# ## old version, which uses a linear KDE for the violin plots (incorrect)
# image <- ggplot() +
#   geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y,
#                                    color = violinwidth), show.legend = FALSE) +
#   # Re-use geom_violin to plot the outline
#   geom_violin(data = rs_boxplot, aes(x = as.integer(QueryPoint), y = Angle, fill = QueryPoint),
#               color = "white", alpha = 0, draw_quantiles = c(0.25, 0.5, 0.75),
#               show.legend = FALSE) +
#   scale_x_continuous(breaks = breaks, labels = labels) +
#   scale_fill_discrete(guide = "none") +
#   scale_color_viridis_c(name="Density") +
#   geom_segment(
#     data = peak_results_df_long,
#     aes(x = as.integer(QueryPoint) - 0.45,
#         xend = as.integer(QueryPoint) + 0.45,
#         y = Angle, yend = Angle, linetype="Peak Angle"),
#     color = "red", linewidth = 1.2, alpha = 0.8
#   ) +
#   geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.Rev, shape = "STICCC"), size = 4, color = "purple", alpha=0.8) +
#   geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.Rev.180, shape = "STICCC"), size = 4, color = "purple", alpha=0.8) +
#   geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.scVelo, shape = "scVelo"), size = 6, color = "#56B4E9", alpha=0.7) +
#   geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.scVelo.Dyn, shape = "scVelo-dyn"), size = 6, color = "#F0E442", fill="#F0E442", alpha=0.7) +
#   geom_point(data = hline_df, aes(x = 1:nrow(hline_df), y = Angle.veloVI, shape = "veloVI"), size = 6, color = "#D55E00", alpha=0.7) +
#   guides(
#     shape = guide_legend(override.aes = list(color = c("purple", "#56B4E9", "#F0E442",  "#D55E00"), size = 4, alpha = 0.8, breaks=c("STICCC", "scVelo", "scVelo-dyn","veloVI")))
#   ) +
#   theme_sticcc() +
#   theme(axis.text.x = element_text(angle=90)) +
#   scale_linetype_manual(name = "Curve", values = c("Peak Angle" = "dotted")) +
#   scale_shape_manual(name = "Method", values = c("STICCC" = 16, "scVelo" = 17, "scVelo-dyn"=25, "veloVI"=15), breaks=c("STICCC", "scVelo", "scVelo-dyn","veloVI")) +  # Different point shapes
#   labs(x = "Trajectory Point", y = "Angle")
# image



