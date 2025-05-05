# Load library
library(readr)
library(pheatmap)
library(dplyr)
library(ggpubr)

# Set working directory
getwd()
setwd("C:/Users/gui_l/OneDrive/Mestrado/Artigo 2")

# Set seed to always have the same results
set.seed(123)

# Load data
library(readxl)
dados <- read_excel("picrust2_data.xlsx", sheet = "Planilha1")

# 
variables <- colnames(dados[,-c(1:5)])
rownames(dados) <- dados$SampleID
dados[,-c(1:5)] <- scale(dados[,-c(1:5)])
dados <- as.data.frame(t(dados))

dados_nl <- dados[-c(1:5),]
dados_nl <- lapply(dados_nl, as.numeric)
dados_nl <- as.data.frame(dados_nl)
dados_nl <- scale(dados_nl)
rownames(dados_nl) <- variables

dados_tg <- dados[c(1:5),]
group <- as.matrix(dados_tg[2,])

# Pheatmap plot
p <- pheatmap(dados_nl,
              cluster_rows = T,
              cluster_cols =  F,
              scale = "row",
              cellwidth = 15,
              cellheigh = 15,
              fontsize = 10,
              border_color = "black",
              color = colorRampPalette(c("#073647", "white", "#A33333"))(11))

svg(filename = 'heatmap.svg', width = 16)
p
dev.off()

pheatmap(dados_nl,
         scale = "row")


# Baloon plot
dados <- read_delim("~/dados/Guilherme/heatmap/dados", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
dados_clay <- dados %>% 
  filter(Soil == "Clay")
dados_clay <- dados_clay[,-c(1:5)]
plot_balloon_clay <- as.data.frame(colMeans(dados_clay))

dados_sandy <- dados %>% 
  filter(Soil == "Sandy clay loam")
dados_sandy<- dados_sandy[,-c(1:5)]
plot_balloon_sandy <- as.data.frame(colMeans(dados_sandy))

plot_balloon <- cbind(plot_balloon_clay, plot_balloon_sandy)
balloon <- ggballoonplot(plot_balloon)

svg(filename = "balloon.svg")
balloon
dev.off()
