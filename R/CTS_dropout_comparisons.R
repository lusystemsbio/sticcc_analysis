rm(list=ls())
source("R/HPCFunctions.R")
source("R/flow_analysis_functions.R")


topoName <- "CTS"
radius <- c(0,1)
angleHow <- "avg-"
ctrl <- "ND" # 1 or "ND"
minMag <- 0


## Compare noise level vs angle change
vic_full <- readRDS(file.path(getwd(), topoName, "data", paste0("vic_",topoName,"_inVectors.Rds")))
if(topoName == "repressilator") {
  subset <- colnames(readRDS(file.path(getwd(), topoName, "data", paste0("vic_",topoName,"_BoolODE_dropout_",0.3,"_may152023.Rds"))))  
} else {
  subset <- colnames(readRDS(file.path(getwd(), topoName, "data", paste0("vic_",topoName,"_BoolODE_dropout_",0.3,"_nov132023.Rds"))))  
}

#subset <- gsub("V","",subset)

colnames(vic_full) <- paste0("V",colnames(vic_full))
vic_wt <- vic_full[,subset]

# invert v2
vic_wt@metadata$vectors_in <- vic_wt@metadata$vectors_in * -1
colData(vic_wt)$dX_in <- -1 * colData(vic_wt)$dX_in
colData(vic_wt)$dY_in <- -1 * colData(vic_wt)$dY_in

wt_output_irrev <- getAnglesInRadius(vic_wt, radius, how = angleHow, minMagnitude = minMag)
#wt_output_irrev$SampleID <- rownames(wt_output_irrev)
wt_output_irrev$Condition = "ND"
comb_df <- wt_output_irrev[,c("gridpoint","arclen","vec_angle","vec_angle_shuffled","Condition", "dir", "dir_shuffled")]

#aList <- c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")[1:10]
aList <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
for(i in seq_along(aList)) {
  print(paste0("i = ",i))
  
  a <- aList[i]
  
  # Read in kd simulation results
  if(topoName == "repressilator") {
    vic_fname <- file.path(getwd(), topoName, "data", paste0("vic_",topoName,"_BoolODE_dropout_",a,"_may152023.Rds"))  
  } else if(topoName == "CTS") {
    vic_fname <- file.path(getwd(), topoName, "data", paste0("vic_",topoName,"_BoolODE_dropout_",a,"_nov132023.Rds"))  
  }
  vic_kd <- readRDS(vic_fname)
  # invert v2
  vic_kd@metadata$vectors_in <- vic_kd@metadata$vectors_in * -1
  colData(vic_kd)$dX_in <- -1 * colData(vic_kd)$dX_in
  colData(vic_kd)$dY_in <- -1 * colData(vic_kd)$dY_in
  
  # Compute angles for kd
  kd_output <- getAnglesInRadius(vic_kd, radius, how = angleHow, minMagnitude = minMag)
  #kd_output$SampleID <- rownames(kd_output)
  kd_output$Condition <- a
  
  # Identify shared points
  overlap <- kd_output$gridpoint[which(kd_output$gridpoint %in% comb_df$gridpoint)]
  comb_df <- comb_df[which(comb_df$gridpoint %in% overlap),]
  kd_output <- kd_output[overlap,c("gridpoint","arclen","vec_angle","vec_angle_shuffled","Condition", "dir", "dir_shuffled")]
  
  # Join dataframes
  comb_df <- rbind(comb_df, kd_output)
  
}


comb_df$Condition <- factor(comb_df$Condition, levels=c("ND", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))

## Compute summary stats on angle based on noise
image <- ggplot(data=comb_df) +
  geom_boxplot(mapping=aes(x=Condition, color=Condition,y=vec_angle)) +
  ggtitle("Absolute angles vs noise")
image


## Compute differences from (WT or a=1)?
#diff_df <- comb_df[,c("gridpoint","vec_angle","arclen","Condition")]
#diff_df <- pivot_wider(diff_df, names_from = "Condition", values_from = c("arclen","vec_angle"))

diff_df_long <- comb_df
wt_arcs <- comb_df[which(comb_df$Condition == ctrl),"arclen"]
wt_angles <- comb_df[which(comb_df$Condition == ctrl),"vec_angle"]
for(a in c("ND",aList)) {
  #diff_df[,paste0("arclen_",a)] <- diff_df[,paste0("arclen_",a)] - diff_df$arclen_WT
  #diff_df[,paste0("vec_angle_",a)] <- diff_df[,paste0("vec_angle_",a)] - diff_df$vec_angle_WT
  diff_df_long[which(diff_df_long$Condition == a),"arclen"] <- 
    diff_df_long[which(diff_df_long$Condition == a),"arclen"] - wt_arcs
  
  
  diff_df_long[which(diff_df_long$Condition == a),"vec_angle"] <- 
    sapply(diff_df_long[which(diff_df_long$Condition == a),"vec_angle"] - wt_angles, minAngleCopy)
  
  diff_df_long[which(diff_df_long$Condition == a),"vec_angle_shuffled"] <- 
    sapply(diff_df_long[which(diff_df_long$Condition == a),"vec_angle_shuffled"] - wt_angles, minAngleCopy)
}



## Compute summary stats on angle based on noise
# 
# image <- ggplot(data=diff_df_long) +
#   geom_boxplot(mapping=aes(x=Condition, color=Condition,y=vec_angle)) +
#   ggtitle(paste0("Angle Change vs noise in ",topoName))
# image
# 
# boxplot_fname <- file.path(getwd(),topoName,"manplots",paste0("BoolODE_boxplot_diffs_",angleHow,"_control=",ctrl,".pdf"))
# pdf(file=boxplot_fname, height = 10, width = 10)
# print(image)
# dev.off()





# boxplots with one null dist
nulldf <- diff_df_long[,c("gridpoint","arclen","vec_angle","Condition")]
nulldist <- diff_df_long[which(diff_df_long$Condition=="ND"),c("gridpoint","arclen","vec_angle_shuffled","Condition")]
nulldist$Condition <- "Null"
colnames(nulldist)[3] <- "vec_angle"
nulldf <- rbind(nulldf, nulldist)
nulldf$Condition <- factor(nulldf$Condition, levels=c("Null","ND","0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))



## Compute summary stats on angle based on noise
image <- ggplot(data=nulldf) +
  geom_boxplot(mapping=aes(x=Condition, color=Condition,y=vec_angle)) +
  ggtitle(paste0("Angle change vs noise in ",topoName)) +
  ylab("Angle Change (Radians)") +
  theme(axis.title = element_text(size=30), axis.text = element_text(size=22))
image

boxplot_fname <- file.path(getwd(),topoName,"manplots",paste0("BoolODE_boxplot_diffs_",angleHow,"_control=",ctrl,"_",Sys.Date(),"_withNullDist.pdf"))
pdf(file=boxplot_fname, height = 10, width = 10)
print(image)
dev.off()



## Compute significant differences

hist(nulldf[which(nulldf$Condition=="WT"),"vec_angle"])

hist(nulldf[which(nulldf$Condition=="0.1"),"vec_angle"])
hist(nulldf[which(nulldf$Condition=="0.2"),"vec_angle"])

t.test(nulldf[which(nulldf$Condition=="WT"),"vec_angle"], nulldf[which(nulldf$Condition=="0.1"),"vec_angle"])

t.test(nulldf[which(nulldf$Condition=="0.1"),"vec_angle"], nulldf[which(nulldf$Condition=="0.2"),"vec_angle"])
t.test(nulldf[which(nulldf$Condition=="0.1"),"vec_angle"], nulldf[which(nulldf$Condition=="0.3"),"vec_angle"])




## Look at WT vs 0.0 drop probability results
vic_fname <- file.path(getwd(), topoName, "data", paste0("vic_",topoName,"_BoolODE_dropout_0.0_feb112023.Rds"))  
vic_kd <- readRDS(vic_fname)
wt_expr <- assay(vic_wt, "normcounts")
noDrop_expr <- assay(vic_kd, "normcounts")

test <- noDrop_expr - 0.001
diff <- test - wt_expr
