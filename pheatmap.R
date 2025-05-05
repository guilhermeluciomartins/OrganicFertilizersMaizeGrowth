# Pheatmap para dados de Picrust2

# Load packages
library(pheatmap)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(readxl)


# set and get the directory work
getwd()
#setwd("C:/Users/gui_l/OneDrive/R")
setwd("C:/Users/gui_l/OneDrive/Mestrado/Artigo 2")

# Load data files
data <- read.table(file = "pheatmap.txt", dec = ".", sep = "\t", header = T)


# Transform multiple columns into z-scores using scale()
data_zscore <- as.data.frame(scale(data[6:21]))

# Add columns 1 to 5 from the original "data" to "data_zscore"
merged_data <- cbind(data[1:5], data_zscore)

# View the transformed data.frame
head(merged_data)
str(merged_data)

# Drop columns
merged_data <- merged_data[, -c(1, 3:5)]

# Transform the data into a matrix
matrix_data <- as.matrix(merged_data)


is.na(matrix_data)
is.infinite(matrix_data)

# Basic heatmap
pheatmap(matrix_data, main = "basic heatmap")

# heatmap with cluster break
pheatmap(matrix_data, cutree_rows = 5, cutree_cols = 3, main = "heatmap with clusterized output")



# load dataset
df = read.csv("https://reneshbedre.github.io/assets/posts/heatmap/hm_data.csv", row.names="Gene")
# convert to matrix
df_mat = data.matrix(df)
