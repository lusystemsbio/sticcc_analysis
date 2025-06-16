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
library(biomaRt) # Gene name conversion
library(Seurat) # cell cycle phase assignment & UMAP
library(zoo) # rollapply
library(tidyr) # pivot_longer
## For cell cycle gene set
library(msigdbr)
library(dplyr)
library(stringr)
library(scran) # cell cycle phase assignment (not really used though now)

set.seed(123)
source("R/sticcc_analysis_utilities.R")

#remotes::install_version("dplyr", version = "1.1.3")
#remotes::install_version("dbplyr", version = "2.3.3")


# global params
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
topoName <- "cellcycle_human_u2os_v3"
assignCCPhases <- FALSE

#unsure if these are needed
forcePrep <- FALSE
forcePCA <- FALSE
forceSTICCC <- FALSE
saveNetworkPlot <- FALSE

nSamples <- 10000
pseudocount <- T
numClusters <- 3


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



##### IMPORT EXPERIMENTAL DATA #####
exprMat_unlog <- t(read.csv(file.path(pythonDataDir, paste0("U2OS_counts_norm.csv")), 
                           row.names = 1, header = T))
exprMat_norm <- log2(exprMat_unlog + 1)
metadata <- data.frame(SampleID=colnames(exprMat_norm))
rownames(metadata) <- metadata$Sample

# ##### CELL CYCLE PHASE ASSIGNMENT #####
# if(assignCCPhases) {
#   metadata <- data.frame(SampleID=colnames(exprMat_norm))
#   rownames(metadata) <- metadata$Sample
#   
#   tmp <- SingleCellExperiment(assays=SimpleList(counts=exprMat_unlog, logcounts=exprMat_norm))
#   cc_genes <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
#   assignments <- cyclone(tmp, pairs = cc_genes)
#   
#   # Add to metadata
#   metadata$SCE_Phase <- assignments$phases
#   metadata$SCE_G1.score <- assignments$scores$G1
#   metadata$SCE_S.score <- assignments$scores$S
#   metadata$SCE_G2M.score <- assignments$scores$G2M
#   
#   saveRDS(metadata, file.path(outputDir,"u2os_metadata.Rds"))
# } else {
#   metadata <- readRDS(file.path(outputDir,"u2os_metadata.Rds"))
# }




##### CONVERT GENE NAMES #####
# Connect to Ensembl BioMart
exprMat_renamed_fname <- file.path(outputDir, "exprMat_norm_renamed.Rds")
if(!file.exists(exprMat_renamed_fname)) {
  #ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
  ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  # Remove version numbers from Ensembl IDs if present
  ensembl_ids <- gsub("\\..*$", "", rownames(exprMat_norm))
  
  # Get gene symbols
  annotations <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Make sure to keep only unique mappings
  annotations <- annotations[annotations$hgnc_symbol != "", ]
  annotations <- annotations[!duplicated(annotations$ensembl_gene_id), ]
  
  # Match and rename
  rownames(exprMat_norm) <- annotations$hgnc_symbol[match(ensembl_ids, annotations$ensembl_gene_id)]
  rownames(exprMat_unlog) <- annotations$hgnc_symbol[match(ensembl_ids, annotations$ensembl_gene_id)]
  
  # Optionally, remove rows with NA rownames (unmatched IDs)
  exprMat_norm <- exprMat_norm[!is.na(rownames(exprMat_norm)), ]
  exprMat_unlog <- exprMat_unlog[!is.na(rownames(exprMat_unlog)), ]
  
  # Remove duplicates
  exprMat_norm <- exprMat_norm[!duplicated(rownames(exprMat_norm)),]
  exprMat_unlog <- exprMat_unlog[!duplicated(rownames(exprMat_unlog)),]
  
  saveRDS(exprMat_norm, exprMat_renamed_fname)
  saveRDS(exprMat_unlog, file.path(outputDir, "exprMat_unlog_renamed.Rds"))
  saveRDS(annotations, file.path(outputDir, "gene_name_conversions.Rds"))
} else {
  exprMat_norm <- readRDS(exprMat_renamed_fname)
  exprMat_unlog <- readRDS(file.path(outputDir, "exprMat_unlog_renamed.Rds"))
  annotations <- readRDS(file.path(outputDir, "gene_name_conversions.Rds"))
}


##### SEURAT CELL CYCLE & PCA #####
# Transpose to match Seurat's expected format (genes as rows, cells as columns)
seurat_fname <- file.path(outputDir,"U2OS_seurat.Rds")
if(!file.exists(seurat_fname)) {
  seurat_obj <- CreateSeuratObject(counts = exprMat_unlog, data=exprMat_norm)
  
  # QC
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 7500 & percent.mt < 30)
  
  # Built-in cell cycle gene lists
  s.genes <- Seurat::cc.genes$s.genes
  g2m.genes <- Seurat::cc.genes$g2m.genes
  
  # Make sure gene symbols match the format (e.g., capitalize if needed)
  common_genes <- rownames(seurat_obj)
  s.genes <- intersect(s.genes, common_genes)
  g2m.genes <- intersect(g2m.genes, common_genes)
  
  # Score cell cycle and assign phase
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
  
  # Extract results
  seurat_metadata <- seurat_obj@meta.data[, c("S.Score", "G2M.Score", "Phase")]
  metadata[colnames(seurat_obj),"Seurat_S.Score"] <- seurat_obj@meta.data$S.Score
  metadata[colnames(seurat_obj),"Seurat_G2M.Score"] <- seurat_obj@meta.data$G2M.Score
  metadata[colnames(seurat_obj),"Seurat_Phase"] <- seurat_obj@meta.data$Phase
  metadata[colnames(seurat_obj),"nCount_RNA"] <- seurat_obj@meta.data$nCount_RNA
  metadata[colnames(seurat_obj),"nFeature_RNA"] <- seurat_obj@meta.data$nFeature_RNA
  
  # 1. Find highly variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 2. Scale data (only variable features)
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  
  # 3. Run PCA/UMAP/TSNE
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- RunUMAP(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- RunTSNE(seurat_obj)
  
  # # Filter outliers?
  # pca_coords <- Embeddings(seurat_obj, "pca")[, 1:10]
  # # Compute Euclidean distances from the median
  # center <- apply(pca_coords, 2, median)
  # distances <- apply(pca_coords, 1, function(x) sqrt(sum((x - center)^2)))
  # # Identify outliers
  # threshold <- median(distances) + 5 * mad(distances)
  # outliers <- which(distances > threshold)
  # # 6. Optionally visualize
  # hist(distances, breaks = 50, main = "Distance from PCA center")
  # abline(v = threshold, col = "red")
  # seurat_obj_filtered <- subset(seurat_obj, cells = setdiff(Cells(seurat_obj), Cells(seurat_obj)[outliers]))
  # 
  saveRDS(seurat_obj, seurat_fname)
  saveRDS(metadata, file.path(outputDir,"u2os_metadata.Rds"))
} else {
  seurat_obj <- readRDS(seurat_fname)
  metadata <- readRDS(file.path(outputDir,"u2os_metadata.Rds"))
}

# Plot PCA
DimPlot(seurat_obj, reduction = "pca", group.by = "Phase")
DimPlot(seurat_obj, reduction = "umap", group.by = "Phase")
DimPlot(seurat_obj, reduction = "tsne", group.by = "Phase")


##### CHECK CELL CYCLE GENE EXPRESSION #####
genes <- unique(c(topo$Source, topo$Target, "CDKN1A", "MCM2","H2AX"))
gene_phases <- list("CCND1" = "G1", # expressed all but separate node
                    "E2F1" = "G1", # expressed in S
                    "E2F7" = "S", # expressed all but G1, noisy
                    "CCNA2" = "G2M", # expressed most in G2M, moderate in S
                    "FZR1" = "", # expressed in G2M
                    "CCNB1" = "G2M", # expressed in G2M
                    "CDC20" = "G2M", # expressed in G2M
                    "CDKN1B" = "", # low, but highest in G2M
                    "CCNE1" = "S", # expressed in S
                    "FBXO5" = "", # expressed S, G2M
                    "RB1" = "") # noisy, inconclusive

tsne_embeddings <- Embeddings(seurat_obj, reduction = "tsne")
tsne_df <- cbind(tsne_embeddings[,c(1:2)], metadata[colnames(seurat_obj),], t(exprMat_norm[genes,colnames(seurat_obj)]))

ggplot(data=tsne_df, aes(x=tSNE_1, y=tSNE_2, shape=Seurat_Phase, color=Seurat_Phase)) +
  geom_point()

##### PCA & STICCC SETUP #####
# create SCE object

reduction.use <- "UMAP"
reduction.use.dist <- "PCA"

stic_fname <- file.path(outputDir,paste0("stic_",topoName,"_dist=",reduction.use.dist,"_plot=",reduction.use,".Rds"))
if(!file.exists(stic_fname) | forcePrep) {
  stic <- sticSE(topo = topo, normData = exprMat_norm[,colnames(seurat_obj)],
                 topoName = topoName, expName = expName)
  
  # add metadata
  #stic <- prepMetadata(stic, exprMat_norm, cluster = F, k = numClusters)
  colData(stic) <- DataFrame(metadata[colnames(seurat_obj),])
  colnames(stic) <- colData(stic)$SampleID
  
  # add PCA
  #reducedDim(stic, "PCA") <- pca_df
  #stic@metadata$pca_data <- list(sdev=NA, rotation=pca_loadings, center=NA, scale=NA)
  #stic@metadata$pca_summary <- NA#summary(pca)
  
  #stic <- runPCA(stic, save=T, overwrite=forcePCA, fname=file.path(outputDir,"PCA_res.Rds"))
  
  # Compute variance and proportion explained
  sdev_seurat <- seurat_obj[["pca"]]@stdev
  pca_var_seurat <- sdev_seurat^2
  pca_var_explained_seurat <- pca_var_seurat / sum(pca_var_seurat)
  
  
  # Add dimension reductions from Seurat
  pca_coords <- Embeddings(seurat_obj, reduction = "pca")
  reducedDim(stic, "PCA") <- pca_coords
  stic@metadata$pca_data <- list(sdev=seurat_obj[["pca"]]@stdev, 
                                rotation=Loadings(seurat_obj, reduction="pca"),
                                center=NA, 
                                scale=F)#pca[1:4]
  stic@metadata$pca_summary <- list(sdev=seurat_obj[["pca"]]@stdev,
                                   rotation=NA,
                                   center=NA,
                                   scale=F,
                                   x=NA,
                                   importance=data.frame(PC1=c(sdev_seurat[1], pca_var_explained_seurat[1]),
                                                         PC2=c(sdev_seurat[2], pca_var_explained_seurat[2])))#summary(pca)
  
  UMAP_coords <- Embeddings(seurat_obj, reduction = "umap")
  reducedDim(stic, "UMAP") <- UMAP_coords
  
  tsne_coords <- Embeddings(seurat_obj, reduction = "tsne")
  reducedDim(stic, "TSNE") <- tsne_coords
  
  #ggplot(reducedDim(stic, "PCA"), aes(x=PC1, y=PC2, color=seurat_metadata$Phase)) +
  #  geom_point(size=3) +
  #  theme_sticcc()
  
  # compute grid based on PCA
  stic@metadata$params$plotDim <- reduction.use
  stic <- computeGrid(stic)
  
  # compute pairwise distance between points
  if(reduction.use.dist == "TSNE" | reduction.use.dist == "UMAP") {
    stic@metadata$params$nDistPCs <- 2
  } else {
    stic@metadata$params$nDistPCs <- 10
  }
  stic <- computeDist(stic, reduction = reduction.use.dist)
  
  saveRDS(stic, stic_fname)
} else {
  stic <- readRDS(stic_fname)
}






# or TSNE
# stic@metadata$params$plotDim <- "TSNE"



##### RUN STICCC #####
# compute trajectories
stic@metadata$topo <- topo
stic@metadata$params$nPCs <- 2


if(!file.exists(stic_fname) | forceSTICCC) {
  #undebug(runSTICCC)
  stic@metadata$params$sample_radius <- 0.2
  stic <- runSTICCC(stic, v2=T, invertV2=T)
  saveRDS(stic, stic_fname)
} else {
  stic <- readRDS(stic_fname)
}
stic <- computeGrid(stic, grid.length = 18)

##### PLOT STICCC VECTORS #####

## First, plot each with individual vectors on PCA
st_vectors_v1 <- stic@metadata$vectors
st_vectors_v2 <- stic@metadata$vectors_in
st_net <- as.data.frame(st_vectors_v1 + st_vectors_v2)
st_rev <- as.data.frame((st_vectors_v1 - st_vectors_v2) / 2)
scalingFactor <- 8

# not actually pca since we're using umap/tsne here
pca_df <- reducedDim(stic,reduction.use)
pca_df <- merge(pca_df, metadata, by="row.names")
rownames(pca_df) <- pca_df$Row.names
pca_df <- pca_df[,-1]
colnames(pca_df)[1:2] <- c("X","Y")
plot_xlab <- paste0(reduction.use,"1")
plot_ylab <- paste0(reduction.use,"2")
xMin <- floor(min(pca_df[,1]))
xMax <- ceiling(max(pca_df[,1]))
yMin <- floor(min(pca_df[,2]))
yMax <- ceiling(max(pca_df[,2]))

net_df <- merge(pca_df, st_net, by="row.names")
image <- ggplot(net_df, aes(x=X, y=Y)) +
  geom_point(aes(color=Seurat_Phase), alpha=0.8, size=4) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("STICCC Net Flow on ",expName," cells")) +
  guides(alpha="none", size="none", color=guide_legend(title = "Phase", override.aes = list(size = 5))) +
  #theme(axis.text = element_text(size=16), axis.title = element_text(size=20))  +
  geom_segment(aes(xend=X+dPC_1*scalingFactor, yend=Y+dPC_2*scalingFactor),
               arrow = arrow(length = unit(0.2,"cm"))) +
  theme_bw() +
  theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) +
  scale_color_manual(values=c(cbPalette))
image

pdf(file = file.path(plotDir, paste0(reduction.use,"_netflow_rad=",stic@metadata$params$sample_radius,".pdf")), width = 10, height = 10)
print(image)
dev.off()



rev_df <- merge(pca_df, st_rev, by="row.names")
image <- ggplot(rev_df, aes(x=X, y=Y)) +
  geom_point(aes(color=Seurat_Phase), alpha=0.8, size=4) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  scale_size(range=c(1.75, 3)) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  ggtitle(paste0("STICCC Reversibility on ",expName," cells")) +
  guides(alpha="none", size="none", color=guide_legend(title = "Phase", override.aes = list(size = 5))) +
  #theme(axis.text = element_text(size=16), axis.title = element_text(size=20))  +
  geom_segment(aes(xend=X+dPC_1*scalingFactor, yend=Y+dPC_2*scalingFactor),
               arrow = arrow(length = unit(0.2,"cm"))) +
  theme_bw() +
  theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) +
  scale_color_manual(values=c(cbPalette))
image
pdf(file = file.path(plotDir, paste0(reduction.use,"_rev_rad=",stic@metadata$params$sample_radius,".pdf")), width = 10, height = 10)
print(image)
dev.off()


stic <- computeGridVectors(stic, unitVectors = F, combine = T, how = "net")
grid.df <- as.matrix(stic@metadata$grid.df)
grid.df[is.nan(grid.df)] <- 0
grid.df <- cbind(grid.df,log(grid.df[,7]+0.1))
colnames(grid.df)[8] <- "logMagnitude"

plot_df <- pca_df

image <- ggplot() +
  geom_point(data=plot_df[,],mapping=aes(x=X,y=Y, color=Seurat_Phase),size=2.4) + 
  theme_bw() +
  metR::geom_streamline(data = grid.df, 
                        aes(x = x.points, y = y.points, dx = dx, dy = dy), 
                        L = 4, res = 20,  arrow.angle=30, n=20, jitter = 4) + 
  guides(alpha="none", size="none", color=guide_legend(title = "Phase", override.aes = list(size = 5))) +
  labs(fill="Density") + 
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) +
  scale_color_manual(values=c(cbPalette))
image  

pdf(file = file.path(plotDir, paste0(reduction.use,"_netflow_rad=",stic@metadata$params$sample_radius,"_streamline.pdf")), width = 10, height = 10)
print(image)
dev.off()


stic <- computeGridVectors(stic, unitVectors = T, combine = T, how = "net")
plotGrid(stic, colorVar="Seurat_Phase", outputDir = plotDir, return = T, minMagnitude = 0.001)


stic <- computeGridVectors(stic, unitVectors = F, combine = T, how = "rev")
grid.df <- as.matrix(stic@metadata$grid.df)
grid.df[is.nan(grid.df)] <- 0
grid.df <- cbind(grid.df,log(grid.df[,7]+0.1))
colnames(grid.df)[8] <- "logMagnitude"

plot_df <- pca_df

image <- ggplot() +
  geom_point(data=plot_df[,],mapping=aes(x=X,y=Y, color=Seurat_Phase),size=2.4) + 
  theme_bw() +
  metR::geom_streamline(data = grid.df, 
                        aes(x = x.points, y = y.points, dx = dx, dy = dy), 
                        L = 4, res = 20,  arrow.angle=30, n=20, jitter = 4) + 
  guides(alpha="none", size="none", color=guide_legend(title = "Phase", override.aes = list(size = 5))) +
  labs(fill="Density") + 
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  xlim(xMin,xMax) +
  ylim(yMin,yMax) +
  theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) +
  scale_color_manual(values=c(cbPalette))
image  

pdf(file = file.path(plotDir, paste0(reduction.use,"_rev_rad=",stic@metadata$params$sample_radius,"_streamline.pdf")), width = 10, height = 10)
print(image)
dev.off()


stic <- computeGridVectors(stic, unitVectors = T, combine = T, how = "rev")
plotGrid(stic, colorVar="Seurat_Phase", outputDir = plotDir, return = T, minMagnitude = 0.001, plotSuffix = "rev")


##### VECTORS BY ANGLE #####

angle_df <- as.data.frame(colData(stic)[,c("X","Y")])
angle_df$Angle <- atan2(angle_df$Y, angle_df$X)

# Adjust angles to be from 0 to 2pi
angle_df$Angle <- (angle_df$Angle + 2*pi) %% (2*pi)

#topoGenes <- unique(c(topo$Source, topo$Target))
traj_df <- as.data.frame(colData(stic)[,c("SampleID", "Seurat_Phase", 
                                          "X", "Y", 
                                          "dX", "dY", "dX_in", "dY_in")])
traj_df$Angle <- angle_df$Angle

traj_df$dX_Net <- (traj_df$dX + traj_df$dX_in) / 2
traj_df$dY_Net <- (traj_df$dY + traj_df$dY_in) / 2

traj_df$dX_Rev <- (traj_df$dX - traj_df$dX_in) / 2
traj_df$dY_Rev <- (traj_df$dY - traj_df$dY_in) / 2

traj_df$Mag_Net <- sqrt(traj_df$dX_Net^2 + traj_df$dY_Net^2)
traj_df$Mag_Rev <- sqrt(traj_df$dX_Rev^2 + traj_df$dY_Rev^2)


traj_df <- traj_df %>%
  rowwise() %>%
  mutate(
    r_norm = sqrt(X^2 + Y^2),
    r_hat_x = X / r_norm,
    r_hat_y = Y / r_norm,
    
    # Orthogonal unit vector for tangential direction
    t_hat_x = -r_hat_y,
    t_hat_y =  r_hat_x,
    
    # Signed radial components
    mag_rad_net = dX_Net * r_hat_x + dY_Net * r_hat_y,
    mag_rad_rev = dX_Rev * r_hat_x + dY_Rev * r_hat_y,
    
    # Signed tangential components (dot with perpendicular vector)
    mag_tan_net = dX_Net * t_hat_x + dY_Net * t_hat_y,
    mag_tan_rev = dX_Rev * t_hat_x + dY_Rev * t_hat_y
  ) %>%
  ungroup()

traj_df$Ratio_Tan <- traj_df$mag_tan_net / traj_df$mag_tan_rev
traj_df$Ratio <- traj_df$Mag_Net / traj_df$Mag_Rev


# Plot with polar coordinates
window_size <- 120
traj_df_states_onehot <- traj_df
traj_df_states_onehot$PhaseS <- ifelse(traj_df_states_onehot$Seurat_Phase == "S", 1, 0)
traj_df_states_onehot$PhaseG1 <- ifelse(traj_df_states_onehot$Seurat_Phase == "G1", 1, 0)
traj_df_states_onehot$PhaseG2M <- ifelse(traj_df_states_onehot$Seurat_Phase == "G2M", 1, 0)

traj_df_ma_states <- traj_df_states_onehot %>%
  arrange(Angle) %>%
  mutate(
    MA_G1 = rollapply(PhaseG1, width = window_size, FUN = mean, partial = TRUE, align = 'center'),
    MA_S = rollapply(PhaseS, width = window_size, FUN = mean, partial = TRUE, align = 'center'),
    MA_G2M = rollapply(PhaseG2M, width = window_size, FUN = mean, partial = TRUE, align = 'center'),
    MA_Tan_Ratio = rollapply(Ratio_Tan, width = window_size, FUN = mean, partial = TRUE, align = 'center',na.rm=T),
    MA_Ratio = rollapply(Ratio, width = window_size, FUN = mean, partial = TRUE, align = 'center',na.rm=T)
  )



# Pivot data to long format for plotting
traj_df_ma_states_long <- traj_df_ma_states[,c("Angle", "MA_G1", "MA_S", "MA_G2M", 
                                               "MA_Tan_Ratio",
                                               "MA_Ratio")] %>%
  pivot_longer(
    cols = -Angle, 
    names_to = "Phase",  # This is the new column that will store the previous column names
    values_to = "Moving_Average"  # This is the new column that will store the median values
  )


# Plotting the moving averages
ggplot(traj_df_ma_states_long, aes(x = Angle, y=Moving_Average)) +
  geom_line(aes(color = Phase), size=2) +
  labs(title = "Dominant Cell Cycle Phase",
       x = "Angle (degrees)",
       y = "Moving Average",
       color = "Legend") +
  theme_sticcc()



# Plotting the moving averages
combined_df <- traj_df_ma_states_long
combined_df[which(combined_df$Phase == "MA_Tan_Ratio"),"Moving_Average"] <- abs(combined_df[which(combined_df$Phase == "MA_Tan_Ratio"),"Moving_Average"])

# scale tangential ratio
max_tan_ratio <- max(combined_df[which(combined_df$Phase == "MA_Tan_Ratio"),"Moving_Average"])
combined_df[which(combined_df$Phase == "MA_Tan_Ratio"),"Moving_Average"] <- combined_df[which(combined_df$Phase == "MA_Tan_Ratio"),"Moving_Average"] / max_tan_ratio

max_ratio <- max(combined_df[which(combined_df$Phase == "MA_Ratio"),"Moving_Average"])
combined_df[which(combined_df$Phase == "MA_Ratio"),"Moving_Average"] <- combined_df[which(combined_df$Phase == "MA_Ratio"),"Moving_Average"] / max_ratio

combined_df[which(combined_df$Phase == "MA_G1"),"Phase"] <- "G1"
combined_df[which(combined_df$Phase == "MA_G2M"),"Phase"] <- "G2M"
combined_df[which(combined_df$Phase == "MA_S"),"Phase"] <- "S"
combined_df[which(combined_df$Phase == "MA_Tan_Ratio"),"Phase"] <- "Net/Rev (Tan, Scaled)"
combined_df[which(combined_df$Phase == "MA_Ratio"),"Phase"] <- "Net/Rev (Scaled)"

combined_df$Phase <- factor(combined_df$Phase, levels = c("G1","S","G2M","Net/Rev (Scaled)", "Net/Rev (Tan, Scaled)"))

# plot full vector MA
image <- ggplot(combined_df[which(combined_df$Phase != "Net/Rev (Tan, Scaled)"),], 
                aes(x = Angle, y = Moving_Average)) +
  geom_line(aes(color = Phase), size = 2) +
  labs(x = "Angle (radians)",
       y = "Proportion of Cells (moving average)",
       color = "Series") +
  scale_x_continuous(
    breaks = seq(0, 2 * pi, by = pi / 2),
    labels = c("0", expression(pi/2), expression(pi),
               expression(3*pi/2), expression(2*pi))
  ) +
  theme_sticcc()

image


pdf(file.path(plotDir,"phases_vs_ratio.pdf"),width = 10, height = 10)
print(image)
dev.off()


# plot tangential vector MA
image <- ggplot(combined_df[which(combined_df$Phase != "Net/Rev (Scaled)"),], 
                aes(x = Angle, y = Moving_Average)) +
  geom_line(aes(color = Phase), size = 2) +
  labs(x = "Angle (radians)",
       y = "Proportion of Cells (moving average)",
       color = "Series") +
  scale_x_continuous(
    breaks = seq(0, 2 * pi, by = pi / 2),
    labels = c("0", expression(pi/2), expression(pi),
               expression(3*pi/2), expression(2*pi))
  ) +
  theme_sticcc()

image


pdf(file.path(plotDir,"phases_vs_ratio_tan.pdf"),width = 10, height = 10)
print(image)
dev.off()



# plot cell-wise ratios
my_breaks <- round(unlist(lapply(c(-2, 0, 2), exp)), 2)
image <- ggplot(traj_df, aes(x=X, y=Y, color=Ratio)) +
  geom_point(size=3) +
  theme_sticcc() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_color_gradient2(name = "t", trans = "log",
                        breaks = my_breaks, labels = my_breaks) +
  labs(#title = "Net Flow/Magnitude Ratio",
       color = "Ratio")
image

pdf(file.path(plotDir,"cellwise_ratio.pdf"),width = 10, height = 10)
print(image)
dev.off()

# plot cell-wise ratios (tan)
my_breaks <- round(unlist(lapply(c(-2, 0, 2), exp)), 2)
traj_df$Ratio_Tan_Abs <- abs(traj_df$Ratio_Tan)
image <- ggplot(traj_df, aes(x=X, y=Y, color=Ratio_Tan_Abs)) +
  geom_point(size=3) +
  theme_sticcc() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_color_gradient2(name = "t", trans = "log",
                        breaks = my_breaks, labels = my_breaks) +
  labs(#title = "Net Flow/Magnitude Ratio (Tan)",
       color = "Ratio")
image

pdf(file.path(plotDir,"cellwise_ratio_tan.pdf"),width = 10, height = 10)
print(image)
dev.off()


