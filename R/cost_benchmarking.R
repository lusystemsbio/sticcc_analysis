##### GLOBALS & SETUP #####
rm(list=ls())
library(STICCC)  # vector inference
library(sRACIPE) # simulate circuits
library(ggplot2) # visualize data
library(microbenchmark) # benchmarking time 
library(dplyr) # dataframe management

topo_name_list <- c("CTS", "repressilator", "cellcycle_v2", "emt_net")
numClusters <- 2 # constant value used here, clustering not relevant to benchmarking results
numTrials <- 3


##### BENCHMARKING #####
# Estimate computational cost (space & time) as a function of 1) num cells and 2) num edges
# 1k, 2k, 3k, 5k, 10k for REP (3-gene, 3-edge), CTS (4-gene, 6-edge), cellcycle sims (), EMT sims ()
sizes <- c(1000, 2000, 3000, 5000, 10000)
results_all <- list()
results_dir <- file.path(getwd(),"benchmarking_summary")
if(!dir.exists(results_dir)) {
  dir.create(results_dir)
}
plot_dir <- file.path(results_dir,"plots")
if(!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
results_all_fname <- file.path(results_dir, "benchmarking_results.Rds")

# Loop over input topos
for(topoName in topo_name_list) {
  # directory setup
  topoDir <- file.path(getwd(),topoName)
  outputDir = file.path(topoDir,"data_benchmarking")

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
  
  
  # Allocate results dataframe
  nEdges <- nrow(topo)
  genes <- unique(c(topo$Source, topo$Target))
  nGenes <- length(genes)
  topo_results <- data.frame(nModels=rep(sizes, each=numTrials), 
                             Trial=rep(1:numTrials, times=length(sizes)), 
                             nEdges=nEdges, 
                             nGenes=nGenes,
                             CPUTime=NA, 
                             Mem=NA,
                             TopoName=topoName
                             )
  topo_results_fname <- file.path(outputDir, paste0(topoName,"_benchmark_results.Rds"))
  
  # Simulate 10k models
  sim_fname <- file.path(outputDir, paste0("simData_",topoName,".Rds"))
  if(!file.exists(sim_fname)) {
    racipe <- simTopo(topo, numModels = max(sizes))
    saveRDS(racipe, sim_fname)
  } else {
    racipe <- readRDS(sim_fname)
  }
  
  # normalize data using built in function, or manually
  exprMat <- assay(racipe)[,]
  racipe_norm <- sracipeNormalize(racipe)
  exprMat_norm <- assay(racipe_norm)[,]
  
  
  # Iterate over sizes
  for(size in sizes) {
    # Downsample
    subset <- sample(1:max(sizes), size)
    exprMat_ds <- exprMat[,subset]
    exprMat_norm_ds <- exprMat_norm[,subset]
    
    # Create STICCC object
    stic <- sticSE(topo = topo, exprMat = exprMat_ds, normData = exprMat_norm_ds,
                   topoName = topoName, expName = paste0(topoName, "_benchmarking_",size))
    # add metadata
    stic <- prepMetadata(stic, exprMat_norm_ds, cluster = T, k = numClusters)
    
    # run PCA
    stic <- runPCA(stic, save=T, overwrite=T, fname = file.path(outputDir,paste0("PCA_res_",size,".Rds")))
    
    # compute grid based on PCA
    stic <- computeGrid(stic)
    
    # compute pairwise distance between points
    stic <- computeDist(stic)
    
    # Run STICCC with profiling
    size_list <- c()
    benchmark_data <- microbenchmark({
      stic <- runSTICCC(stic, v2=T, invertV2=T)
      size_list <- c(size_list, object.size(stic))
    }, times = numTrials)
    
    # Store data
    for(i in 1:numTrials) {
      topo_results[which(topo_results$nModels == size & topo_results$Trial == i),"CPUTime"] <- benchmark_data$time[i]
      topo_results[which(topo_results$nModels == size & topo_results$Trial == i),"Mem"] <- object.size(stic)
    }
    
    saveRDS(topo_results, topo_results_fname)
    
    
  }
  
  results_all[[topoName]] <- topo_results
  saveRDS(results_all, results_all_fname)
}




##### PLOT RESULTS #####
benchmark_data <- do.call('rbind',readRDS(results_all_fname))

## Dataset size vs runtime
# Summarize data: calculate mean and standard error of CPUTime
summary_data <- benchmark_data %>%
  group_by(TopoName, nModels) %>%
  summarise(
    mean_cpu = mean(CPUTime),
    se_cpu = sd(CPUTime) / sqrt(n()),
    mean_mem = mean(Mem),
    se_mem = sd(Mem) / sqrt(n()),
    .groups = "drop"
  )

summary_data$mean_cpu <- summary_data$mean_cpu / (1e9)
summary_data$se_cpu <- summary_data$se_cpu / (1e9)
summary_data[which(summary_data$TopoName == "cellcycle_v2"),"TopoName"] <- "Cell Cycle"
summary_data[which(summary_data$TopoName == "repressilator"),"TopoName"] <- "REP"
summary_data[which(summary_data$TopoName == "emt_net"),"TopoName"] <- "EMT"

image <- ggplot(summary_data, aes(x = nModels, y = mean_cpu, color = TopoName, group = TopoName)) +
  geom_line(linewidth=1.5) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = mean_cpu - se_cpu, ymax = mean_cpu + se_cpu), width = 0.2) +
  labs(
    x = "Number of Models",
    y = "CPU Time (s, mean ± SE)"#,
    #title = "CPU Time vs Number of Models by Topology"
  ) +
  theme_minimal() +
  theme(axis.text = element_text(size=28), axis.title = element_text(size=36))
image


pdf(file.path(plot_dir,"Time_vs_Models.pdf"), width = 10, height = 10)
print(image)
dev.off()

## Dataset size vs memory
image <- ggplot(summary_data, aes(x = nModels, y = mean_mem, color = TopoName, group = TopoName)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_mem - se_mem, ymax = mean_mem + se_mem), width = 0.2) +
  labs(
    x = "Number of Models",
    y = "Memory Usage (bytes, mean ± SE)",
    title = "Memory Usage vs Number of Models by Topology"
  ) +
  theme_minimal()
image


pdf(file.path(plot_dir,"Memory_vs_Models.pdf"), width = 10, height = 10)
print(image)
dev.off()

## Number of edges vs runtime
summary_data_2 <- benchmark_data %>%
  group_by(nEdges, nModels) %>%
  summarise(
    mean_cpu = mean(CPUTime),
    se_cpu = sd(CPUTime) / sqrt(n()),
    mean_mem = mean(Mem),
    se_mem = sd(Mem) / sqrt(n()),
    .groups = "drop"
  )
summary_data_2$nModels <- factor(summary_data_2$nModels)

image <- ggplot(summary_data_2, aes(x = nEdges, y = mean_cpu, color = nModels, group = nModels)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_cpu - se_cpu, ymax = mean_cpu + se_cpu), width = 0.2) +
  labs(
    x = "Number of Edges",
    y = "CPU Time (ns, mean ± SE)",
    title = "CPU Time vs Number of Edges"
  ) +
  theme_minimal()
image


pdf(file.path(plot_dir,"Time_vs_Edges.pdf"), width = 10, height = 10)
print(image)
dev.off()

## Number of edges vs memory
image <- ggplot(summary_data_2, aes(x = nEdges, y = mean_mem, color = nModels, group = nModels)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_mem - se_mem, ymax = mean_mem + se_mem), width = 0.2) +
  labs(
    x = "Number of Edges",
    y = "Memory (bytes, mean ± SE)",
    title = "Memory vs Number of Edges"
  ) +
  theme_minimal()
image


pdf(file.path(plot_dir,"Memory_vs_Edges.pdf"), width = 10, height = 10)
print(image)
dev.off()

## Number of nodes vs runtime
summary_data_3 <- benchmark_data %>%
  group_by(nGenes, nModels) %>%
  summarise(
    mean_cpu = mean(CPUTime),
    se_cpu = sd(CPUTime) / sqrt(n()),
    mean_mem = mean(Mem),
    se_mem = sd(Mem) / sqrt(n()),
    .groups = "drop"
  )
summary_data_3$nModels <- factor(summary_data_3$nModels)
#summary_data_3$nGenes <- factor(summary_data_3$nGenes)
summary_data_3$mean_cpu <- summary_data_3$mean_cpu / (1e9)
summary_data_3$se_cpu <- summary_data_3$se_cpu / (1e9)

image <- ggplot(summary_data_3, aes(x = nGenes, y = mean_cpu, color = nModels, group = nModels)) +
  geom_line(linewidth=1.5) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = mean_cpu - se_cpu, ymax = mean_cpu + se_cpu), width = 0.2) +
  labs(
    x = "Number of Genes",
    y = "CPU Time (s, mean ± SE)"#,
    #title = "CPU Time vs Number of Genes"
  ) +
  theme_minimal() +
  theme(axis.text = element_text(size=28), axis.title = element_text(size=36))
image


pdf(file.path(plot_dir,"Time_vs_Genes.pdf"), width = 10, height = 10)
print(image)
dev.off()

## Number of nodes vs memory
image <- ggplot(summary_data_3, aes(x = nGenes, y = mean_mem, color = nModels, group = nModels)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_mem - se_mem, ymax = mean_mem + se_mem), width = 0.2) +
  labs(
    x = "Number of Genes",
    y = "Memory (bytes, mean ± SE)",
    title = "Memory vs Number of Genes"
  ) +
  theme_minimal() +
  theme(axis.text = element_text(size=28), axis.title = element_text(size=36))
image


pdf(file.path(plot_dir,"Memory_vs_Genes.pdf"), width = 10, height = 10)
print(image)
dev.off()












