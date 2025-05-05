### data analysis with microeco package - agricultural waste paper - Guilherme Martins (16/07/2022)

# load the package
library(microeco)

# Get and set directory
getwd()
setwd("C:/Users/GuilhermeL/OneDrive/Mestrado/Artigo 2")
path <- "C:/Users/GuilhermeL/OneDrive/Mestrado/Artigo 2"
list.files(path)



# load the dataset
# metadata table
sample_info_16S <- read.table("sample_info_KO.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# feature table
otu_table_16S <- read.table("otu_table_KO2.txt", header = TRUE, row.names = 1, sep = "\t", dec = ",")

# taxonomic assignment table
taxonomy_table_16S <- read.table("taxonomy_table_KO2.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# load the environmental data table if it is not in sample table
env_data_16S <- read.table("env_data_KO.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# use pipe operator in magrittr package
library(magrittr) # needs to be run every time you start R and want to use %>%

# set.seed is used to fix the random number generation to make the results repeatable
set.seed(123)

# make the plotting background same with the tutorial
library(ggplot2)
theme_set(theme_bw())

# Make sure that the data types of sample_table, otu_table and tax_table are all data.frame format as the following part shows
class(otu_table_16S)
otu_table_16S[1:5, 1:5]

class(taxonomy_table_16S)
taxonomy_table_16S[1:5, 1:1]





#######################################################################
            ### Organize data for further analysis ###


# Remove NAs to proceed with diversity analysis
# make the taxonomic information unified, very important
taxonomy_table_16S %<>% tidy_taxonomy

# make sure that the rownames of sample information table are sample names
class(sample_info_16S)
sample_info_16S[1:5, ]

# make sure that the environmental data are stored as env data
class(env_data_16S)
env_data_16S[1:5, 1:5]

# In R6 class, '$new' is the original method used to create a new object of class
# Let's create a microtable object with more information
dataset <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S)
dataset



################################################################################

# Define color palette
c1 <- c("gray30", "tomato2", "tomato4", "seagreen2", "seagreen4", "steelblue2", "steelblue4", "gray70")

# Reorder groups
dataset$sample_table$Group %<>% factor(., levels = c("CT", "FC", "PL", "CM"))
dataset$sample_table$Type %<>% factor(., levels = c("Unfertilized", "Fertilized", "Fresh", "Composted"))
str(dataset$sample_table)


################################################################################
### Composition-based class ###


# taxa_abund list in the object of microtable class must be first calculated
# create trans_abund object



# use 10 Phyla with the highest abundance in the dataset.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Funct", ntaxa = 13)

# Relative abundance with two facets at Phylum level
# require package ggh4x, first run install.packages("ggh4x") if not installed
#install.packages("ggh4x")
t1$plot_bar(others_color = "grey70",
            facet = c("Group", "Type"),
            xtext_keep = FALSE,
            legend_text_italic = FALSE,
            barwidth = 1
            )

# save plot with 600 dpi resolution
dev.print(tiff, "relative_abundance_functions.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# show the heatmap with the high abundant genera
# show 30 taxa at Genus level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Funct", ntaxa = 48)

# Plot
t1$plot_heatmap(facet = c("Group", "Type"),
                xtext_keep = FALSE,
                withmargin = TRUE,
                plot_colorscale = "log10",
                color_values = rev(RColorBrewer::brewer.pal(n = 9, name = "Blues"))
                )

#  save plot with 600 dpi resolution
dev.print(tiff, "heatmap_functions.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")




# The trans_venn class is used for venn analysis, i.e. shared and unique taxa
# To analyze the unique and shared OTUs of groups, we first merge samples according to the "Group" column of sample_table

# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Type")
# dataset1 is a new microtable object

# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")

# Plot
t1$plot_venn() # The integer is OTU number # The percentage data is the sequence number/total sequence number

# Save plot with 600 dpi resolution
dev.print(tiff, "venn_diagram.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# transform the results of venn plot to the traditional feature-sample table, that is, another object of microtable class
dataset1 <- dataset$merge_samples(use_group = "Type")
t1 <- trans_venn$new(dataset1)





##########################################################################
###  Model-based class ###

# the differential test result is stored in the object$res_diff
# metastat analysis at Genus level
t1 <- trans_diff$new(dataset = dataset,
                     method = "metastat",
                     group = "Type",
                     taxa_level = "Lv3"
                     )

# t1$res_diff is the differential test result
# t1$res_abund is the group abundance

# the metastat can run the comparisons for each paired group
# the user should use select_group to select the required pair


# select_group should be one of groups in t1$res_diff$Comparison
# FC_f vs FC_c differential abundance
t1$plot_diff_abund(use_number = 1:20,
                   col = c("tomato2", "tomato4"),
                   select_group = "Mineral - Composted"
                   )

#  save plot with 600 dpi resolution
dev.print(tiff, "FC_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# LEfSe combines the non-parametric test and linear discriminant analysis
t1 <- trans_diff$new(dataset = dataset,
                     method = "lefse",
                     group = "Type",
                     alpha = 0.01,
                     lefse_subgroup = NULL
                     )

# see t1$res_diff for the result

# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 1)

#  save plot with 600 dpi resolution
dev.print(tiff, "LEfSE_functions.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# ANOVA method and transpose
#install.packages("agricolae")
t1 <- trans_diff$new(dataset = dataset,
                     method = "anova",
                     group = "Type",
                     taxa_level = "Funct",
                     filter_thres = 0.001
                     )

# Plot
t1$plot_diff_abund(use_number = 1:25,
                   col = c1,
                   add_sig = T,
                   coord_flip = F
                   )

#  save plot with 600 dpi resolution
dev.print(tiff, "ANOVA_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=35, units="cm")


# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset,
                     method = "KW",
                     group = "Type",
                     taxa_level = "Lv3",
                     filter_thres = 0.001
                     )

# Plot
t1$plot_diff_abund(use_number = 1:15,
                   col = c1,
                   add_sig = T,
                   coord_flip = F
                   )

#  save plot with 600 dpi resolution
dev.print(tiff, "KW_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")






#################################################################################
### Explainable class ###

# If there may be some NA in the user's env data,add complete_na = TRUE when creating the trans_env object

# add_data is used to add the environmental data
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)

# show the autocorrelations among variables

# use group parameter to show the distributions of variables and the autocorrelations across groups
#install.packages("GGally")
t1$cal_autocor(group = "Group", col = c1)

#  save plot with 600 dpi resolution
dev.print(tiff, "autocorr.tiff", compression = "lzw", res=600, height=40, width=40, units="cm")


# RDA at Genus level
t1$cal_ordination(method = "RDA", taxa_level = "Funct")
t1$res_ordination
t1$res_ordination_R2
t1$res_ordination_terms
t1$res_ordination_axis

#r.squared adj.r.squared 
#0.8400305     0.6933918 

#Df   Variance       F Pr(>F)    
#Group     1 0.00035530  5.7094  0.001 ***
#Type      1 0.00029710  4.7742  0.003 ** 
#Unified   1 0.00055181  8.8672  0.001 ***
#N         1 0.00068687 11.0375  0.001 ***
#P         1 0.00040582  6.5212  0.001 ***
#K         1 0.00028336  4.5534  0.009 ** 
#Ca        1 0.00018139  2.9149  0.025 *  
#Mg        1 0.00006122  0.9837  0.383    
#S         1 0.00051319  8.2465  0.001 ***
#Cu        1 0.00001868  0.3001  0.668    
#Zn        1 0.00056669  9.1063  0.001 ***
#Residual 12 0.00074677                   
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.

t1$trans_ordination(show_taxa = 10,
                    adjust_arrow_length = TRUE,
                    max_perc_env = 1.5,
                    max_perc_tax = 1.5,
                    min_perc_env = 0.2,
                    min_perc_tax = 0.2
                    )

# t1$res_rda_trans is the transformed result for plotting

# Plot
t1$plot_ordination(plot_color = "Group",
                   plot_shape = "Type",
                   col = c1,
                   plot_type = c("point")
                   )

#  save plot with 600 dpi resolution
dev.print(tiff, "RDA_genus.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# Mantel test can be used to check significant correlations between environmental variables and distance matrix
t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
t1$res_mantel


# perform a correlation heatmap between environmental and taxa data at Phylum level
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
t1$cal_cor(use_data = "Lv2", cor_method = "spearman", p_adjust_method = "fdr")
# return t1$res_cor

# plot the correlation results using plot_cor function
# filter phylum that do not have at least one ***
t1$plot_cor()
#  save plot with 600 dpi resolution
dev.print(tiff, "functions_correlation2.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# Pheatmap method
# clustering heatmap; require pheatmap package
# Let's take another color pallete: 
#install.packages('pheatmap')

# Plot
t1$plot_cor(pheatmap = TRUE,
            color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
            )

#  save plot with 600 dpi resolution
dev.print(tiff, "pheatmap_phylum_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



#correlations between environmental variables and taxa for different groups
# calculate correlations for different groups using parameter by_group

# first create trans_diff object as a demonstration
#install.packages("randomForest")
t2 <- trans_diff$new(dataset = dataset,
                     method = "rf",
                     group = "Type",
                     rf_taxa_level = "Lv2"
                     )

# then create trans_env object
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)

# use other_taxa to select taxa you need
t1$cal_cor(by_group = "Type",
           use_data = "other",
           cor_method = "spearman", 
           p_adjust_method = "fdr",
           other_taxa = t2$res_diff$Taxa[1:25]
           )

# Plot
t1$plot_cor()

#  save plot with 600 dpi resolution
dev.print(tiff, "env_corr.tiff", compression = "lzw", res=600, height=15, width=40, units="cm")
