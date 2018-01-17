# Import data
# setting up working directory
setwd("C:\\MalignancyLearning")
getwd()
full_data <- read.csv("data/data_wisconsin.csv")
# importing caret
# install.packages("caret")
require(caret)
# importing ggplot2 library for visualization
library(ggplot2)
library(reshape2)
library(mlbench)
library(corrplot)
library(pcaGoPromoter)
library(ellipse)
library(gridExtra)
library(grid)
# sneak peek of data
head(full_data)
# number of rows in dataset
nrow(full_data)
# Breakdown for every column
summary(full_data)
# Quick debreif of dataframe ('str' = structure)
str(full_data)
# random variables distributed uniformly 
# runif() 

# Identify features that most likely cause tumors to be malignant 
# applying filter will allow us to assess all the cases that were deemed malignant  
mal_filter <- full_data$diagnosis == 'M'
malignancy_features <- full_data[mal_filter,]
# Making sure it worked
head(malignancy_features, n = 15)
nrow(malignancy_features)

# Alright let's plot some information to gain some understanding how
# all those parameters are associated to each other
# Var1 <- colnames(malignancy_features)
# Var2 <- rownames(malignancy_features)
plot1 <- ggplot(data=full_data, aes(x=smoothness_mean, y=area_mean, colour = diagnosis))
plot1 + geom_point() +
  ggtitle("Relation b/w Area and smoothness of cells") +
  xlab("Area") +
  ylab("Smoothness") +
  theme(axis.title.x = element_text(colour = "Red", size=20),
        axis.title.y = element_text(colour = "Blue", size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=10),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        plot.title = element_text(colour="Black",size=20, family=""))
# 
pca_func <- function(data, groups, title, print_ellipse = TRUE) {
  
  # perform pca and extract scores
  pcaOutput <- pca(data, printDropped = FALSE, scale = TRUE, center = TRUE)
  pcaOutput2 <- as.data.frame(pcaOutput$scores)
  
  # define groups for plotting
  pcaOutput2$groups <- groups
  
  # when plotting samples calculate ellipses for plotting (when plotting features, there are no replicates)
  if (print_ellipse) {
    
    centroids <- aggregate(cbind(PC1, PC2) ~ groups, pcaOutput2, mean)
    conf.rgn  <- do.call(rbind, lapply(unique(pcaOutput2$groups), function(t)
      data.frame(groups = as.character(t),
                 ellipse(cov(pcaOutput2[pcaOutput2$groups == t, 1:2]),
                         centre = as.matrix(centroids[centroids$groups == t, 2:3]),
                         level = 0.95),
                 stringsAsFactors = FALSE)))
    
    plot <- ggplot(data = pcaOutput2, aes(x = PC1, y = PC2, group = groups, color = groups)) + 
      geom_polygon(data = conf.rgn, aes(fill = groups), alpha = 0.2) +
      geom_point(size = 2, alpha = 0.6) + 
      scale_color_brewer(palette = "Set1") +
      labs(title = title,
           color = "",
           fill = "",
           x = paste0("PC1: ", round(pcaOutput$pov[1], digits = 2) * 100, "% variance"),
           y = paste0("PC2: ", round(pcaOutput$pov[2], digits = 2) * 100, "% variance"))
    
  } else {
    
    # if there are fewer than 10 groups (e.g. the predictor classes) I want to have colors from RColorBrewer
    if (length(unique(pcaOutput2$groups)) <= 10) {
      
      plot <- ggplot(data = pcaOutput2, aes(x = PC1, y = PC2, group = groups, color = groups)) + 
        geom_point(size = 2, alpha = 0.6) + 
        scale_color_brewer(palette = "Set1") +
        labs(title = title,
             color = "",
             fill = "",
             x = paste0("PC1: ", round(pcaOutput$pov[1], digits = 2) * 100, "% variance"),
             y = paste0("PC2: ", round(pcaOutput$pov[2], digits = 2) * 100, "% variance"))
      
    } else {
      
      # otherwise use the default rainbow colors
      plot <- ggplot(data = pcaOutput2, aes(x = PC1, y = PC2, group = groups, color = groups)) + 
        geom_point(size = 2, alpha = 0.6) + 
        labs(title = title,
             color = "",
             fill = "",
             x = paste0("PC1: ", round(pcaOutput$pov[1], digits = 2) * 100, "% variance"),
             y = paste0("PC2: ", round(pcaOutput$pov[2], digits = 2) * 100, "% variance"))
      
    }
  }
  
  return(plot)
  
}


p1 <- pca_func(data = t(full_data[, 3:32]), groups = as.character(bc_data_2$diagnosis), title = "Breast cancer dataset 2: Samples")
p2 <- pca_func(data = full_data[, 3:32], groups = as.character(colnames(bc_data_2[, 3:32])), title = "Breast cancer dataset 2: Features", print_ellipse = FALSE)
grid.arrange(p1, p2, ncol = 2, widths = c(0.4, 0.6))



# Correlation matrix to find attributes that are highly correlated, and rank features based on importance
# set.seed(7)

corr_matrix <- cor(full_data[,3:32])
corrplot(corr_matrix,order='hclust')
# From the correlation plot, we can see some of the obvious correlations, radius strongly correlates with area, perimeter.
# 
highCorr <- findCorrelation(corr_matrix, cutoff = 0.5,verbose=T)
colnames(full_data[,highCorr])
# highCorr
str(malignancy_features)
mal_high_corr <- cor(malignancy_features[,3:32])
highCorrMal <- findCorrelation(mal_high_corr,cutoff=0.7,verbose=T)
corrplot(mal_high_corr)

# Ranking feature importance
control <- trainControl(method="repeatedcv", number = 10, repeats = 10)
imp_features <- function(model, title){
  # estimating importance
  importance <- varImp(model, scale=TRUE)

  # prepare dataframes for plotting
  importance_df_1 <- importance$importance
  importance_df_1$group <- rownames(importance_df_1)
  importance_df_2 <- importance_df_1
  importance_df_2$Overall <- 0

  importance_df <- rbing(importance_df_1, importance_df_2)

  plot <- ggplot() +
    geom_point(data = importance_df_1, aes(x=Overall, y = group, color = group), size = 2) +
    geom_path(data = importance_df, aes(x=Overall, y= group, color = group, group = group),size=1) +
    labs(
      x = "Importance",
      y = "",
      title = title,
      subtitle = "Scaled feature importance",
      caption = "\nDetermined with Random Forest and repeated cross validation (repeats-10, 10 times)"
    )
  return(plot)
}



