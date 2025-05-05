# Heatmap for PICRUSt2 analysis - Guilherme Martins (14-9-2023)

# Load packages
library(tidyverse)
library(reshape2)
library(vegan)
library(corrplot)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(compositions)
library(readxl)

# set and get the directory work
getwd()
#setwd("C:/Users/gui_l/OneDrive/R")
setwd("C:/Users/GuilhermeL/OneDrive/R")

# Load data files
data <- read.table(file = "heatmap4.txt", dec = ".", sep = "\t", header = T)

# View the transformed data.frame
head(data)
str(data)

# Data for the balloon plot
db_n <- read_excel("balloon_plot.xlsx", sheet = "nitrogen")
db_m <- read_excel("balloon_plot.xlsx", sheet = "methane")
db_c <- read_excel("balloon_plot.xlsx", sheet = "carbon")



# Data transformation ----------------------------------------------------------

# Select the columns to be transformed
clr_data <- data[, 6:21]
# Perform CLR transformation
clr_data <- clr(clr_data)
# Return columns from the original "data"
clr_data <- cbind(data[1:5], clr_data)


# Select rows for the nitrogen, methane and carbon cycles
df_nitrogen <- clr_data[c(1:27), ]
df_methane <- clr_data[c(28:83), ]
df_carbon <- clr_data[c(84:153), ]


# Transform multiple columns into z-scores using scale()
data_zscore_n <- as.data.frame(scale(df_nitrogen[6:21]))
data_zscore_m <- as.data.frame(scale(df_methane[6:21]))
data_zscore_c <- as.data.frame(scale(df_carbon[6:21]))


# Add columns 1 to 5 from the original "data" to "data_zscore"
merged_data_n <- cbind(df_nitrogen[1:5], data_zscore_n)
merged_data_m <- cbind(df_methane[1:5], data_zscore_m)
merged_data_c <- cbind(df_carbon[1:5], data_zscore_c)


# Separate the data frames for each soil type
df_nitrogen_1 <- merged_data_n[, c(1:13)]
df_nitrogen_2 <- merged_data_n[, c(1:5, 14:21)]

df_methane_1 <- df_methane[, c(1:13)]
df_methane_2 <- df_methane[, c(1:5, 14:21)]

df_carbon_1 <- df_carbon[, c(1:13)]
df_carbon_2 <- df_carbon[, c(1:5, 14:21)]


# setting a theme for your plots -----------------------------------------------

tema1 <- theme_bw() +
  theme(
    text = element_text(colour = "gray25"),
    plot.title = element_text(size = 14), 
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = NA),
    axis.text.y.left = element_text(size = 12),
    #axis.text.y.left = element_blank(),
    axis.text.x = element_text(angle = 90, size = 12, vjust = 0.5, hjust = 0.5),
    #axis.ticks.x = element_blank(),
    legend.key.size = unit(0.7,"cm"),
    legend.title = element_text(size = 12)
  )

tema2 <- theme_bw() +
  theme(
    text = element_text(colour = "gray25"),
    plot.title = element_text(size = 14), 
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = NA),
    #axis.text.y.left = element_text(size = 12),
    axis.text.y.left = element_blank(),
    axis.text.x = element_text(angle = 90, size = 12, vjust = 0.5, hjust = 0.5),
    #axis.ticks.x = element_blank(),
    legend.key.size = unit(0.7,"cm"),
    legend.title = element_text(size = 12),
  )



# Reshape the data from wide to long format ------------------------------------

# Nitrogen metabolism ----------------------------------------------------------
melted_nitrogen_1 <- melt(df_nitrogen_1,
                        id.vars = c("Metabolism", "Function", "Gene", "Detail"),
                        measure.vars = c("CT_unfert1","CT_Fert1",
                                         "FC_fresh1", "FC_composted1",
                                         "PL_fresh1", "PL_composted1",
                                         "CM_fresh1", "CM_composted1"))


melted_nitrogen_2 <- melt(df_nitrogen_2,
                          id.vars = c("Metabolism", "Function", "Gene", "Detail"),
                          measure.vars = c("CT_unfert2","CT_Fert2",
                                           "FC_fresh2", "FC_composted2",
                                           "PL_fresh2", "PL_composted2",
                                           "CM_fresh2", "CM_composted2"))


# Plot your graph (clay soil)
p1 = ggplot(melted_nitrogen_1, aes(x = variable, y = Function)) + 
  geom_tile(aes(fill = value, height = 0.95, width = 0.95)) + 
  labs(x = '', y = 'Gene abundance (clr transformed)') +
  scale_fill_gradient2(low="red", mid = "white" , high="blue")  + 
  ggtitle("Clay soil") +
  tema1
p1


# Plot your graph (Sandy Clay Loam)
p11 = ggplot(melted_nitrogen_2, aes(x = variable , y = Function)) + 
  geom_tile(aes(fill = value, height = 0.95, width = 0.95)) + 
  labs(x = '', y = '') +
  scale_fill_gradient2(low="red", mid = "white" , high="blue")  + 
  ggtitle("Sandy clay loam soil") +
  tema2
p11


# Balloon plot
b1 <- ggballoonplot(db_n, x= "Axis", y= "Function" , size = "Mean_clay")
b2 <- ggballoonplot(db_n, x= "Axis", y= "Function" , size = "Mean_sandy")+
  theme(axis.text.y.left = element_blank())


# Merge heatmap plots
nitrogen_figure <- ggarrange(p1, p11, 
                    labels = c("A", "B"),
                    common.legend = T,
                    legend = "right",
                    widths = c(2,1),
                    heights = c(1,1),
                    ncol = 2, nrow = 1)

annotate_figure(nitrogen_figure,
                fig.lab = "Nitrogen metabolism",
                fig.lab.face = "bold",
                fig.lab.size = 14,
                fig.lab.pos = "top.left")

# Saved as SVG (w = 1450 x h = 650)


# Merge balloon plots
balao1 <- ggarrange(b1, b2, 
                             common.legend = F,
                             legend = "right",
                             ncol = 2, nrow = 1)

balao1 # Saved as SVG (w = 1450 x h = 650)



# Methane metabolism -----------------------------------------------------------
melted_methane_1 <- melt(df_methane_1,
                          id.vars = c("Metabolism", "Function", "Gene", "Detail"),
                          measure.vars = c("CT_unfert1","CT_Fert1",
                                           "FC_fresh1", "FC_composted1",
                                           "PL_fresh1", "PL_composted1",
                                           "CM_fresh1", "CM_composted1"))


melted_methane_2 <- melt(df_methane_2,
                          id.vars = c("Metabolism", "Function", "Gene", "Detail"),
                          measure.vars = c("CT_unfert2","CT_Fert2",
                                           "FC_fresh2", "FC_composted2",
                                           "PL_fresh2", "PL_composted2",
                                           "CM_fresh2", "CM_composted2"))


# Plot your graph (clay soil)
p2 = ggplot(melted_methane_1, aes(x = variable , y = Function)) + 
  geom_tile(aes(fill = value, height = 0.95, width = 0.95)) + 
  labs(x = '', y = 'Gene abundance (clr transformed)') +
  #scale_fill_gradient(low="white", high="#317ec2ff")  +
  scale_fill_gradient2(low="red", mid = "white" , high="blue")  + 
  ggtitle("Clay soil") +
  tema1
p2 #+ coord_flip()


# Plot your graph (Sandy clay loam soil)
p22 = ggplot(melted_methane_2, aes(x = variable , y = Function)) + 
  geom_tile(aes(fill = value, height = 0.95, width = 0.95)) + 
  labs(x = '', y = '') +
  #scale_fill_gradient(low="white", high="#317ec2ff")  +
  scale_fill_gradient2(low="red", mid = "white" , high="blue")  + 
  ggtitle("Sandy clay loam soil") +
  tema2
p22 #+ coord_flip()


# Balloon plot
b3 <- ggballoonplot(db_m, x= "Axis", y= "Function" , size = "Mean_clay")
b4 <- ggballoonplot(db_m, x= "Axis", y= "Function" , size = "Mean_sandy")+
  theme(axis.text.y.left = element_blank())


# Merge heatmap plots
methane_figure <- ggarrange(p2, p22, 
                             labels = c("A", "B"),
                             common.legend = T,
                             legend = "right",
                             widths = c(1.8,1),
                             heights = c(1,1),
                             ncol = 2, nrow = 1)

annotate_figure(methane_figure,
                fig.lab = "Methane metabolism",
                fig.lab.face = "bold",
                fig.lab.size = 14,
                fig.lab.pos = "top.left")

# saved as SVG (1450 x 550) (w x h)


# Merge balloon plots
balao2 <- ggarrange(b3, b4, 
                    common.legend = F,
                    legend = "right",
                    ncol = 2, nrow = 1)

balao2 # Saved as SVG (w = 1450 x h = 550)



# Carbon fixation metabolism ----------------------------------------------------
melted_carbon_1 <- melt(df_carbon_1,
                         id.vars = c("Metabolism", "Function", "Gene", "Detail"),
                         measure.vars = c("CT_unfert1","CT_Fert1",
                                          "FC_fresh1", "FC_composted1",
                                          "PL_fresh1", "PL_composted1",
                                          "CM_fresh1", "CM_composted1"))


melted_carbon_2 <- melt(df_carbon_2,
                         id.vars = c("Metabolism", "Function", "Gene", "Detail"),
                         measure.vars = c("CT_unfert2","CT_Fert2",
                                          "FC_fresh2", "FC_composted2",
                                          "PL_fresh2", "PL_composted2",
                                          "CM_fresh2", "CM_composted2"))


# Plot your graph (clay soil)
p3 = ggplot(melted_carbon_1, aes(x = variable , y = Function)) + 
  geom_tile(aes(fill = value, height = 0.95, width = 0.95)) + 
  labs(x = '', y = 'Gene abundance (z-score transformed)') +
  scale_fill_gradient2(low="red", mid = "white" , high="blue")  + 
  ggtitle("Clay soil") +
  tema1
p3


# Plot your graph (Sandy clay loam soil)
p33 = ggplot(melted_carbon_2, aes(x = variable , y = Function)) + 
  geom_tile(aes(fill = value, height = 0.95, width = 0.95)) + 
  labs(x = '', y = '') +
  scale_fill_gradient2(low="red", mid = "white" , high="blue")  + 
  ggtitle("Sandy clay loam soil") +
  tema2
p33


# Balloon plots
b5 <- ggballoonplot(db_c, x= "Axis", y= "Function" , size = "Mean_clay")
b6 <- ggballoonplot(db_c, x= "Axis", y= "Function" , size = "Mean_sandy")+
  theme(axis.text.y.left = element_blank())


# Merge heatmap plots
carbon_figure <- ggarrange(p3, p33, 
                            labels = c("A", "B"),
                            common.legend = T,
                            legend = "right",
                            widths = c(1.4,1),
                            heights = c(1,1),
                            ncol = 2, nrow = 1)

annotate_figure(carbon_figure,
                fig.lab = "Carbon fixation metabolism",
                fig.lab.face = "bold",
                fig.lab.size = 14,
                fig.lab.pos = "top.left")


# saved as SVG (1400 x 450) (w x h)


# Merge balloon plots
balao3 <- ggarrange(b5, b6, 
                    common.legend = T,
                    legend = "right",
                    ncol = 2, nrow = 1)

balao3 # Saved as SVG (w = 1450 x h = 550)


# Merge all plots in one figure ------------------------------------------------

ggarrange(p3, p2, p1, 
          labels = c("A", "B", "C"),
          common.legend = F,
          legend = "right",
          ncol = 1, nrow = 3)

