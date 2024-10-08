

theme_sticcc <- function() {
  font <- "Helvetica"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 28,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 16),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = ,                 #font size
        hjust = 16),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 28),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 22),                #font size
      
      legend.title=element_text(size=18), 
      legend.text=element_text(size=16),
      
      
      
      #axis.text.x = element_text(            #margin for axis text
      #  margin=margin(5, b = 10))
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}


getTrajectoryNeighbors <- function(trajectory,
                                   max_dist,
                                   queryData,
                                   neighborhoodRadius) {
  # Find neighborhood members
  lim_dist <- max_dist * neighborhoodRadius
  est_k <- nrow(trajectory)*neighborhoodRadius*2
  neighbors <- FNN::get.knnx(data=trajectory[,which(colnames(trajectory) != "Time")], query=queryData, k=est_k)
  neighbor_idx <- neighbors$nn.index[which(neighbors$nn.dist <= lim_dist)]
  
  return(neighbor_idx)
  
}



# Calculate Euclidean distance between two vectors
calculate_distance <- function(vec1, vec2) {
  # Calculate the sum of squared differences
  sum_squared_differences <- sum((vec1 - vec2)^2)
  
  # Return the square root of the sum of squared differences
  sqrt(sum_squared_differences)
}

# Find the average angle difference across all points in queryTrajectory, for a set of lags and trajectory
v_obs_along_path <- function(trajectory, # PCA plus Time column
                             lag, # numeric - time gap, single value OR vector the length of queryTrajectory (computed numerically, so should correspond to times, not indices!)
                             sce,
                             queryTrajectory, # either rownames of trajectory, or a dataframe of same ncol to be compared to it
                             neighborhoodRadius,
                             v_pred) {
  
  # Get query data
  if(all(queryTrajectory %in% rownames(trajectory))) {
    querySet <- trajectory[queryPoint,which(colnames(trajectory) != "Time")]  
  } else if(ncol(queryTrajectory) != (ncol(trajectory)-1)) {
    print("Error: queryPoint should be a rowname of trajectory or a vector of length ncol(trajectory)-1")
    
  } else {
    querySet <- queryTrajectory
  }
  
  if(length(lag) == 1) {
    lag <- rep(lag, nrow(querySet))
  } else if (length(lag) != nrow(querySet)) {
    print("Error: lag should be a numeric of length 1 or the length of the queryTrajectory")
    return(NULL)
  }
  
  
  # Prep output data
  summary_data <- data.frame(QueryPoint = rownames(querySet),
                             X = querySet[,1],
                             Y = querySet[,2],
                             Lag=lag,
                             Pred.Angle = NA,
                             NumNeighbors = NA,
                             dX.Obs.Mean=NA,
                             dX.Obs.Median=NA,
                             dY.Obs.Mean=NA,
                             dY.Obs.Median=NA,
                             Obs.Angle.1 = NA,
                             Obs.Mag.1 = NA,
                             Obs.Angle.2 = NA,
                             Obs.Mag.2 = NA,
                             Num.Peaks = NA
  )
  
  # out_df will hold data for violin plots
  out_df <- data.frame()
  
  
  # Iterate over points in querySet
  for(i in seq_along(1:nrow(querySet))) {
    if(i %% 10 == 0) {
      print(paste0("Computing fit for query point ", i, " of ", nrow(querySet)))  
    }
    
    queryData <- querySet[i,]
    
    
    # Get initial points (neighbors of queryPoint)
    neighbors <- getTrajectoryNeighbors(trajectory = trajectory,
                                        max_dist = sce@metadata$max_dist,
                                        queryData = queryData,
                                        neighborhoodRadius = neighborhoodRadius)
    
    # check that there are sufficient neighbors
    if(length(neighbors) < 5) {
      angle_diff <- NA
      next
    }
    
    
    # Leave padding around simulation start and end corresponding to max(lags)
    rownames(trajectory) <- round(trajectory$Time, 1)
    init_idx_padded <- neighbors[which(trajectory[neighbors,"Time"] > lag[i] & trajectory[neighbors,"Time"] < (max(trajectory$Time)-lag[i]))]
    init_states <- trajectory[init_idx_padded,]
    times <- init_states$Time
    
    
    
    # Get final states
    final_times <- round((times + lag[i]), 1)
    final_states <- trajectory[as.character(final_times),]
    median_final_state <- colMedians(as.matrix(final_states), useNames = T)
    mean_final_state <- colMeans(as.matrix(final_states))
    
    # Derive v_obs and convert to angle
    v_obs <- final_states - init_states
    v_obs_median <- median_final_state - queryData
    v_obs_mean <- mean_final_state - queryData
    
    # Check: t should be constant in this new df
    angles <- angle_conversion(as.data.frame(v_obs[,c(1:2)]))
    angle_med <- angle_conversion(as.data.frame(v_obs_median[,c(1:2)]))
    #angle_mean <- angle_conversion(as.data.frame(v_obs_mean[,c(1:2)]))
    
    # Identify peak angles
    angle_hist <- hist(angles, plot=F, breaks=20)
    angle_freqs <- angle_hist$counts
    angle_peaks <- pracma::findpeaks(angle_freqs)
    if(nrow(angle_peaks) == 1) {
      angle_peaks <- as.matrix(t(angle_peaks[order(angle_peaks[,1], decreasing = T),])) # sort in descending order
    } else {
      angle_peaks <- angle_peaks[order(angle_peaks[,1], decreasing = T),] # sort in descending order
    }
    
    angle_peak_vals <- angle_hist$breaks[angle_peaks[,2]]
    angle_peak_mags <- angle_peaks[,1] - min(angle_hist$counts)
    #peak_data[which(peak_data$Lag == lag), "NPeaks"] <- length(angle_peaks)
    #peak_data[which(peak_data$Lag == lag), "Peaks"] <- angle_peaks
    
    
    
    
    
    
    # Store data
    boxplot_data <- data.frame(QueryPoint = rep(rownames(querySet)[i], length(init_idx_padded)),
                               T_Init=times, 
                               T_Final=final_times,
                               Lag=lag[i],
                               Idx=init_idx_padded,
                               Angle=angles
    )
    
    
    
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"Pred.Angle"] <- angle_conversion(v_pred[i,])
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"NumNeighbors"] <- length(init_idx_padded)
    
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"dX.Obs.Mean"] <- v_obs_mean[,1]
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"dX.Obs.Median"] <- v_obs_median[,1]
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"dY.Obs.Mean"] <- v_obs_mean[,2]
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"dY.Obs.Median"] <- v_obs_median[,2]
    
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"Obs.Angle.1"] <- angle_peak_vals[1]
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"Obs.Mag.1"] <- angle_peak_mags[1]
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"Obs.Angle.2"] <- angle_peak_vals[2]
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"Obs.Mag.2"] <- angle_peak_mags[2]
    summary_data[which(summary_data$QueryPoint == rownames(querySet)[i]),"Num.Obs.Peaks"] <- nrow(angle_peaks)
    
    
    
    
    # Add to out_df
    out_df <- rbind(out_df, boxplot_data)
    
    #Pred.Angle = rep(angle_conversion(v_pred[i,]), length(lags)),
    #NumNeighbors = rep(length(init_idx_padded), length(lags)),
    
    
    
    
    
    
    
  }
  
  return(list(Summary=summary_data, BoxplotData=out_df))
  
}




# Modified from orignial function to look for nearest neighbors not within a grid square but by distance
vobs_var_by_t <- function(trajectory, # PCA plus Time column
                          lags, # a vector of time gaps (computed numerically, so should correspond to times, not indices!)
                          sce,
                          queryPoint, # either a rowname of trajectory, or a vector to be compared to it
                          neighborhoodRadius,
                          #smoothingRadius,
                          #grid.df, # from computeGrid
                          #point, # one of the values in grid.df$GridPoint
                          plot=T, # whether to display plot
                          save=T, # whether to save plot
                          fname=NA # filename/path to save plot
                          
) {
  
  # Logic
  if(save & is.na(fname)) {
    print("Warning: No file will be saved because fname is unspecified.")
    save=F
  }
  
  # Get query data
  if(all(queryPoint %in% rownames(trajectory))) {
    queryData <- trajectory[queryPoint,which(colnames(trajectory) != "Time")]  
  } else if(length(queryPoint) != (ncol(trajectory)-1)) {
    print("Error: queryPoint should be a rowname of trajectory or a vector of length ncol(trajectory)-1")
    
  } else {
    queryData <- queryPoint
  }
  
  
  # Get initial points (neighbors of queryPoint)
  neighbors <- getTrajectoryNeighbors(trajectory = trajectory,
                                      max_dist = sce@metadata$max_dist,
                                      queryData = queryData,
                                      neighborhoodRadius = neighborhoodRadius)
  
  # check that there are sufficient neighbors
  if(length(neighbors) < 5) {
    print(paste0("Error: not enough cells within the given neighborhoodRadius, only found ",length(neighbors)))
    return(NULL)
  }
  
  
  # Leave padding around simulation start and end corresponding to max(lags)
  rownames(trajectory) <- round(trajectory$Time, 1)
  init_idx_padded <- neighbors[which(trajectory[neighbors,"Time"] > max(lags) & trajectory[neighbors,"Time"] < (max(trajectory$Time)-max(lags)))]
  init_states <- trajectory[init_idx_padded,]
  times <- init_states$Time
  
  
  # Prep output data
  var_data <- data.frame(Lag = lags,
                         Var = NA,
                         RMSD = NA)
  
  # Iterate over lags
  for(lag in lags) {
    
    # Get final states
    final_times <- round((times + lag), 1)
    final_states <- trajectory[as.character(final_times),]
    
    # Var
    final_state_var <- sum(diag(var(final_states[,which(colnames(final_states) != "Time")])))
    var_data[which(var_data$Lag == lag), "Var"] <- final_state_var
    
    # RMSD
    # Apply the distance calculation for each paired row in df1 and df2
    # Using `mapply` to apply the function across the rows of the two dataframes
    distances <- mapply(calculate_distance, init_states[,which(colnames(init_states) != "Time")], 
                        final_states[,which(colnames(final_states) != "Time")], SIMPLIFY = TRUE)
    
    # Calculate RMSD (Root Mean Squared Distance) across all distances
    rmsd <- sqrt(mean(distances^2))
    var_data[which(var_data$Lag == lag), "RMSD"] <- rmsd
    
  }
  
  
  # Plot variance vs lag
  image <- ggplot(data=var_data) +
    geom_point(aes(x=Lag, y=RMSD), size=3) +
    geom_line(aes(x=Lag, y=RMSD)) +
    ggtitle(paste0("RMSD by time lag, point ",rownames(queryData)[1])) +
    theme_sticcc()
  
  
  
  if(plot) {
    print(image)
  }
  
  if(save) {
    pdf(paste0(fname,".pdf"))
    print(image)
    dev.off()
  }
  
  return(var_data)
  
  
  
}


trajectorySmoothVectors <- function(trajectory, # trajectory in PCA coordinates
                                    sce,
                                    neighborhoodRadius = 0.05,
                                    invertV2 = T,
                                    vec.use = "net") {
  
  out_df <- data.frame(point=1:nrow(trajectory), x=NA, y=NA, dx=NA, dy=NA)
  
  for(pt in seq_along(1:nrow(trajectory))) {
    
    if(pt %% 10 == 0) {
      print(paste0("Computing smooth vector for point ", pt, " of ", nrow(trajectory)))
    }
    
    ptPCA <- trajectory[pt,]
    
    
    v_smooth <- smoothVector(sce = sce,
                             ptPCA,
                             neighborhoodRadius = neighborhoodRadius,
                             invertV2 = T)
    
    
    if(is.null(v_smooth)) {
      out_df[pt,"x"] <- ptPCA[,1]
      out_df[pt,"y"] <- ptPCA[,2]
      out_df[pt,"dx"] <- NA
      out_df[pt,"dy"] <- NA
      
      next
    }
    
    if(vec.use =="v1") {
      v_smooth_pred <- data.frame(dx=v_smooth$v1_x, dy=v_smooth$v1_y)
    } else if (vec.use == "v2") {
      v_smooth_pred <- data.frame(dx=v_smooth$v2_x, dy=v_smooth$v2_y)
    } else if (vec.use == "net") {
      v_smooth_pred <- data.frame(dx=v_smooth$net_x, dy=v_smooth$net_y)
    } else if (vec.use == "rev") {
      v_smooth_pred <- data.frame(dx=v_smooth$rev_x, dy=v_smooth$rev_y)
    }
    
    out_df[pt,"x"] <- ptPCA[,1]
    out_df[pt,"y"] <- ptPCA[,2]
    out_df[pt,"dx"] <- v_smooth_pred[,1]
    out_df[pt,"dy"] <- v_smooth_pred[,2]
    
  }
  
  return(out_df)
}

# Convert the first two components to angles (in radians)
angle_conversion <- function(v) {
  atan2(v[, 2], v[, 1])
}



interpolate_vectors <- function(v1, v2, n) {
  if (length(v1) != length(v2)) {
    stop("Vectors v1 and v2 must have the same length")
  }
  
  # Create a matrix where each column represents linearly spaced points between
  # the corresponding elements of v1 and v2
  interpolated_matrix <- sapply(1:length(v1), function(i) {
    seq(from = v1[i], to = v2[i], length.out = n)
  })
  
  # Convert the matrix to a dataframe
  interpolated_df <- as.data.frame(interpolated_matrix)
  
  return(interpolated_df)
}


plot_v_vs_trajectory <- function(sce,
                                 queryPoint,
                                 neighborhoodRadius = 0.1,
                                 trajectory,
                                 lag, # can be either a numeric or vector of length up to 3
                                 v_pred,
                                 scale_factor = 0.1
) {
  
  # Logic
  # if(save & is.na(fname)) {
  #   print("Warning: No file will be saved because fname is unspecified.")
  #   save=F
  # }
  
  # Get query data
  if(all(queryPoint %in% rownames(trajectory))) {
    queryData <- trajectory[queryPoint,which(colnames(trajectory) != "Time")]  
  } else if(length(queryPoint) != (ncol(trajectory)-1)) {
    print("Error: queryPoint should be a rowname of trajectory or a vector of length ncol(trajectory)-1")
    
  } else {
    queryData <- queryPoint
  }
  
  
  # Get initial points (neighbors of queryPoint)
  neighbors <- getTrajectoryNeighbors(trajectory = trajectory,
                                      max_dist = sce@metadata$max_dist,
                                      queryData = queryData,
                                      neighborhoodRadius = neighborhoodRadius)
  
  # check that there are sufficient neighbors
  if(length(neighbors) < 5) {
    print(paste0("Error: not enough cells within the given neighborhoodRadius, only found ",length(neighbors)))
    return(NULL)
  }
  
  
  # Leave padding around simulation start and end corresponding to max(lags)
  init_idx_padded <- neighbors[which(trajectory[neighbors,"Time"] > max(lag) & trajectory[neighbors,"Time"] < (max(trajectory$Time)-max(lag)))]
  init_states <- trajectory[init_idx_padded,]
  times <- init_states$Time
  
  
  
  # Get final states
  if(length(lag) == 1) {
    final_times <- round((times + lag), 1)
    final_states <- trajectory[which(round(trajectory$Time, 1) %in% final_times),]
    
    
    # Plot vs ensemble distribution
    image <- ggplot(data=reducedDim(sce,"PCA")) +
      geom_point(alpha=0.4, color="gray", aes(x=PC1, y=PC2)) + # ensemble distribution
      geom_point(data=init_states, aes(x=PC1,y=PC2),color="darkblue") + # initial states
      #geom_point(data=final_states, aes(x=PC1,y=PC2),color="lightblue") + # final states
      geom_density2d(data=final_states, aes(x=PC1,y=PC2), alpha=0.5) + # final states
      geom_segment(aes(x=queryData$PC1, # predicted vector
                       y=queryData$PC2,
                       xend=queryData$PC1+v_pred$dx*scale_factor, 
                       yend=queryData$PC2+v_pred$dy*scale_factor), color="black", arrow = arrow(length = unit(0.4,"cm"))) 
    
  } else {
    #colorPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    final_times <- round((times + lag[1]), 1)
    final_states <- trajectory[which(round(trajectory$Time, 1) %in% final_times),]
    final_states$Lag <- lag[1]
    for(lagT in lag[2:length(lag)]) {
      new_final_times <- round((times + lagT), 1)
      new_final_states <- trajectory[which(round(trajectory$Time, 1) %in% new_final_times),]
      new_final_states$Lag <- lagT
      final_states <- rbind(final_states, new_final_states)
    }
    
    # Plot vs ensemble distribution
    image <- ggplot(data=reducedDim(sce,"PCA")) +
      geom_point(alpha=0.4, color="gray", aes(x=PC1, y=PC2)) + # ensemble distribution
      geom_point(data=init_states, aes(x=PC1,y=PC2),color="black") + # initial states
      geom_point(data=final_states, aes(x=PC1,y=PC2, color=Lag), alpha=0.5) + # final states
      #geom_density2d_filled(data=final_states, aes(x=PC1,y=PC2), alpha=0.5) + # final states
      scale_color_gradient(guide="legend") + 
      geom_segment(aes(x=queryData$PC1, # predicted vector
                       y=queryData$PC2,
                       xend=queryData$PC1+v_pred$dx*scale_factor, 
                       yend=queryData$PC2+v_pred$dy*scale_factor), color="black", arrow = arrow(length = unit(0.4,"cm"))) 
    
  }
  
  
  return(image)
  
  
}


plotObsVPred <- function(sce,
                         model
) {
  
  # Create sign vector
  topo <- sce@metadata$topo
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'), warn_missing = FALSE)
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]
  
  # Compute sampling radius
  sampling_radius <- sce@metadata$params$sample_radius * sce@metadata$max_dist
  
  # get posMat (position matrix - this should be either identical or closely correlated w/ the data used to compute dist_mat)
  if(is.null(sce@metadata$params$nPCs)) {
    nPCs <- nrow(sce)
    sce@metadata$params$nPCs <- nPCs
    
  } else {
    nPCs <- sce@metadata$params$nPCs
  }
  
  # use first nPCs PCs as spatial coordinates unless genes are specified
  pcaMat <- reducedDim(sce, "PCA")[colnames(sce),1:nPCs]
  if(sce@metadata$params$useGenes) {
    posMat <- as.data.frame(t(assay(sce, "normcounts")))[colnames(sce),]
  } else {
    posMat <- pcaMat
  }
  
  
  
  posList <- colnames(posMat)
  rownames(posMat) <- colnames(sce)
  
  # get exprMat (normalized expression matrix - high-dimensional data (rows=features, cols=samples))
  exprMat <- assay(sce,"normcounts")
  
  # Begin velocity calculation
  if(sce@metadata$params$verbose) {
    print("Computing inferred velocity")
  }
  
  sample_list <- colData(sce)$SampleID
  colData(sce)$numNeighbors <- NA
  colData(sce)$selfCor <- NA
  colData(sce)$X <- NA
  colData(sce)$Y <- NA
  colData(sce)$dX <- NA
  colData(sce)$dY <- NA
  colData(sce)$fewNeighbors <- FALSE
  colData(sce)$nonInvertible <- FALSE
  colData(sce)$selfCorNA <- FALSE
  
  weighted_vector_list <- vector(mode = "list", length = length(sample_list))
  
  stepi <- 1
  query_point <- model
  
  
  
  ## Get query point data
  query_data <- posMat[query_point,]
  colData(sce)[query_point,"X"] <- pcaMat[query_point,1]
  colData(sce)[query_point,"Y"] <- pcaMat[query_point,2]
  
  ## Find neighbors
  neighbors <- which(sce@metadata$dist_mat[query_point,colnames(sce)] <= sampling_radius)
  dists <- sce@metadata$dist_mat[query_point,neighbors]
  dists <- dists[which(names(dists) != query_point)]
  #neighbors <- neighbors[which(neighbors %in% colnames(sce))]
  
  ## Skip sample if number of neighbors too small
  if(length(neighbors) <= (sce@metadata$params$minNeighbors+1)) {
    colData(sce)[query_point,"fewNeighbors"] <- TRUE
    colData(sce)[query_point,"numNeighbors"] <- length(neighbors)-1
    next
  }
  
  ## Select neighbors
  subset_models <- as.data.frame(posMat[neighbors,])
  subset_models <- subset_models[-which(rownames(subset_models) %in% query_point),]
  colData(sce)[query_point,"numNeighbors"] <- nrow(subset_models)
  
  ## Multiple regression
  # Create relative positional matrix of neighbor samples
  subset_models[,1:ncol(subset_models)] <- apply(subset_models[,1:ncol(subset_models)], 2, listToNumeric)
  deriv_df <- sweep(subset_models, 2, as.numeric(query_data), "-")
  
  # get relative expression matrix & transpose (neighboring cells only)
  geneex_mat <- as.matrix(deriv_df[,posList])
  geneex_mat_transposed <- t(geneex_mat)
  
  # Multiply relative positional matrix by its transpose
  mat_product <- geneex_mat_transposed %*% geneex_mat
  
  # Take the inverse if possible - if not, skip this sample
  inverse_mat <- tryCatch({
    solve(mat_product)
  }, error = function(e) {
    "non-invertible"
  })
  
  if(inverse_mat[1] == "non-invertible") {
    colData(sce)[query_point,"nonInvertible"] <- TRUE
    next
  }
  
  # Multiply the inverted matrix by the transposed positional matrix
  x_mat <- inverse_mat %*% geneex_mat_transposed
  
  
  ## Calculate delayed correlation
  # Self correlation - skip this sample if NA
  src_activity <- as.numeric(exprMat[topo$Source,query_point] * sign_vector)
  self_target_activity <- as.numeric(exprMat[topo$Target,query_point]) 
  self_cor <- cor(self_target_activity, src_activity)
  if(is.na(self_cor)) {
    colData(sce)[query_point,"selfCorNA"] <- TRUE
    next
  } 
  colData(sce)[query_point,"selfCor"] <- self_cor
  
  # For each neighbor:
  deriv_df$DelayedCorr <- NA
  for(model in rownames(deriv_df)) {
    # Correlate neighbor target activity with query source activity * sign_vector
    neighbor_target_activity <- as.numeric(exprMat[topo$Target, model]) 
    neighbor_cor <- cor(neighbor_target_activity,src_activity)
    deriv_df[model,"DelayedCorr"] <- (neighbor_cor - self_cor)
  }
  
  # Remove NA values caused by zero variance vectors
  na_models <- which(is.na(deriv_df$DelayedCorr))
  if(length(na_models) > 0) {
    deriv_df <- deriv_df[-na_models,]
    subset_models <- subset_models[-na_models,]
    x_mat <- x_mat[,-na_models]
  }
  
  
  # Compute product of earlier matrix product and delayed corr vector
  corr_vec <- as.matrix(deriv_df[,"DelayedCorr"])
  b_vec <- x_mat %*% corr_vec
  
  ## Save vector to master list
  saved_vector <- as.data.frame(t(b_vec))
  colnames(saved_vector) <- paste0("d",posList)
  weighted_vector_list[[stepi]] <- saved_vector
  names(weighted_vector_list)[[stepi]] <- query_point
  
  ## Update colData - if all genes were used, convert to pcs for plotting
  if(sce@metadata$params$useGenes) {
    pca <- sce@metadata$pca_data
    pcaVector <- scale(as.data.frame(t(b_vec)), pca$center, pca$scale) %*% pca$rotation
    colData(sce)[query_point,"dX"] <- pcaVector[1]
    colData(sce)[query_point,"dY"] <- pcaVector[2]
  } else {
    colData(sce)[query_point,"dX"] <- saved_vector[1]
    colData(sce)[query_point,"dY"] <- saved_vector[2]
  }
  
  
  
  # loop through neighbors again to compute regression accuracy
  deriv_df$Err <- NA
  deriv_df$Obs_CCC <- NA
  deriv_df$Pred_CCC <- NA
  deriv_df$Dist <- NA
  for(model in rownames(deriv_df)) {
    
    rel_ccc <- deriv_df[model,"DelayedCorr"]
    deriv_df[model,"Obs_CCC"] <- rel_ccc
    pred_ccc <- as.numeric(deriv_df[model,c(1:length(posList))]) %*% as.numeric(b_vec[1:length(posList),])
    deriv_df[model,"Pred_CCC"] <- pred_ccc
    err <- (abs(rel_ccc - pred_ccc))^2
    deriv_df[model,"Error"] <- err
    
    deriv_df[model,"Dist"] <- dists[model]
  }
  colname <- paste0("Error",sce@metadata$params$sample_radius)
  colData(sce)[query_point,colname] <- median(deriv_df$Error, na.rm=T)
  
  for(gene in posList) {
    colname <- paste0("beta_",gene)
    colData(sce)[query_point,colname] <- b_vec[gene,]
  }
  
  
  #plot(deriv_df$Obs_CCC, deriv_df$Pred_CCC, col=deriv_df$Dist)
  image <- ggplot(data=deriv_df, aes(x=Obs_CCC, y=Pred_CCC)) +
    geom_point(aes(color=Dist), size=4) +
    theme(axis.text = element_text(size=24), axis.title = element_text(size=24),
          title = element_text(size=20)) 
  print(image)
  
  
}

listToNumeric <- function(x) {
  return(as.numeric(unname(unlist(x))))
}



