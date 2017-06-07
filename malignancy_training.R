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
# sneak peek of data
head(full_data)
# number of rows in dataset
nrow(full_data)
summary(full_data)
