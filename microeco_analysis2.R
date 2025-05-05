### data analysis with microeco package - agricultural waste paper - Guilherme Martins (16/07/2022)

# load the package
library(microeco)

# Get and set work directory
getwd()
setwd("C:/Users/glmartins/OneDrive/Mestrado/Artigo 2")
path <- "C:/Users/glmartins/OneDrive/Mestrado/Artigo 2"
list.files(path)

# load the example data; 16S rRNA gene amplicon sequencing dataset
# metadata table
sample_info_16S <- read.table("sample_info_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# feature table
otu_table_16S <- read.table("otu_table_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# taxonomic assignment table
taxonomy_table_16S <- read.table("taxonomy_table_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# load the environmental data table if it is not in sample table
env_data_16S <- read.table("env_data_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

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
taxonomy_table_16S[1:5, 1:3]


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

# remove OTUs which are not assigned in the Kingdom "k__Archaea" or "k__Bacteria"
dataset$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
dataset

# make the OTU and sample information consistent across all files in the dataset object
dataset$tidy_dataset()
print(dataset)

# let's use sample_sums() to check the sequence numbers in each sample
dataset$sample_sums() %>% range

# Rarefy data to reduce bias in diversity index
dataset$rarefy_samples(sample.size = 2781)

# Check the sequence numbers again
dataset$sample_sums() %>% range

# let's calculate the taxa abundance at each taxonomic rank
# use default parameters
dataset$cal_abund()

# return dataset$taxa_abund
class(dataset$taxa_abund)

# show part of the relative abundance at Phylum level
dataset$taxa_abund$Phylum[1:5, 1:5]

# The function save_abund() can be used to save the taxa abundance file to a local place easily.
dataset$save_abund(dirpath = "taxa_abund")



##########################################################################
### Define color palette ###
c1 <- c("gray30","gray70", "tomato2", "tomato4", "seagreen2", "seagreen4", "steelblue2", "steelblue4")
c2 <- c("gray50", "tomato3", "seagreen3", "steelblue3")

# reorder groups
dataset$sample_table$Group %<>% factor(., levels = c("CT", "FC", "PL", "CM"))
dataset$sample_table$Type %<>% factor(., levels = c("Unfert", "Fert", "Fresh", "Composted"))
dataset$sample_table$Unified %<>% factor(., levels = c("Unfert", "Fert", "FC_f", "FC_c", "PL_f", "PL_c", "CM_f", "CM_c"))
str(dataset$sample_table)



#######################################################################
### Calculate alpha and beta diversity ###

### Alpha diversity ###

# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = TRUE)

# return dataset$alpha_diversity
class(dataset$alpha_diversity)

# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")



### Calculate beta diversity ###
# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")

# we first create an trans_beta object
# measure parameter can invoke the distance matrix in dataset$beta_diversity
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")

# use PCoA as an example, PCA or NMDS is also available
t1$cal_ordination(ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
t1$plot_ordination(plot_color = "Group", plot_shape = "Type", plot_type = c("point", "ellipse"), col = c1, shape_values = c(15, 16, 17, 18))


################################################################################
### Composition-based class ###

# taxa_abund list in the object of microtable class must be first calculated
# create trans_abund object


# use 10 Phyla with the highest abundance in the dataset.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Family", ntaxa = 12)
# Relative abundance with two facets at Phylum level
# require package ggh4x, first run install.packages("ggh4x") if not installed
#install.packages("ggh4x")
t1$plot_bar(others_color = "grey70", facet = "Group", facet2 = "Type", xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)
# save plot with 600 dpi resolution
dev.print(tiff, "relative_abundance_facet2.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# show the heatmap with the high abundant genera
# show 30 taxa at Genus level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Family", ntaxa = 15)
t1$plot_heatmap(facet = "Group", xtext_keep = F, withmargin = T)
#  save plot with 600 dpi resolution
dev.print(tiff, "heatmap_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# The trans_venn class is used for venn analysis, i.e. shared and unique taxa
# To analyze the unique and shared OTUs of groups, we first merge samples according to the "Group" column of sample_table

# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Group")
# dataset1 is a new microtable object
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn(col = c2) # The integer is OTU number # The percentage data is the sequence number/total sequence number
#  save plot with 600 dpi resolution
dev.print(tiff, "venn_diagram.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# transform the results of venn plot to the traditional feature-sample table, that is, another object of microtable class
dataset1 <- dataset$merge_samples(use_group = "Type")
t1 <- trans_venn$new(dataset1)

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)


# Venn diagram by group
dataset1 <- dataset$merge_samples(use_group = "Group")
# dataset1 is a new microtable object
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn() # The integer is OTU number # The percentage data is the sequence number/total sequence number
#  save plot with 600 dpi resolution
dev.print(tiff, "venn_diagram2.tiff", compression = "lzw", res=600, height=15, width=15, units="cm")


################################################################################
### Diversity-based class ###

# The data_alpha is used for the following differential test and plotting
t1 <- trans_alpha$new(dataset = dataset, group = "Unified")
# return t1$data_stat
t1$data_stat[1:5, ]

# test the differences among groups using:
# Kruskal-Wallis Rank Sum Test (overall test when groups > 2)
t1$cal_diff(method = "KW")
# return t1$res_diff
t1$res_diff[1:5, ]

# Dunn's Kruskal-Wallis Multiple Comparisons (for paired groups when groups > 2)
t1$cal_diff(method = "KW_dunn")
# return t1$res_diff
t1$res_diff[1:5, ]

# ANOVA with multiple comparisons
t1$cal_diff(method = "anova")
# return t1$res_diff
t1$res_diff[1:5, ]

# plot the mean and se of alpha diversity for each group, and add the anova result
t1$cal_diff(method = "anova")
t1$plot_alpha(measure = "Chao1")
t1$plot_alpha(measure = "ACE")
t1$plot_alpha(measure = "Simpson")
t1$plot_alpha(measure = "Shannon")
t1$plot_alpha(measure = "Observed")


# Clustering plot is also a frequently used method

# plot and compare the group distances
# calculate and plot sample distances within groups
t1$cal_group_distance()
# return t1$res_group_distance
t1$plot_group_distance(distance_pair_stat = TRUE)

# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$plot_group_distance(distance_pair_stat = TRUE)

# use replace_name to set the label name, group parameter used to set the color
t1$plot_clustering(group = "Unified", replace_name = c("Group", "Type"))
#  save plot with 600 dpi resolution
dev.print(tiff, "clustering.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# PERMANOVA is often used in the differential test of distances among groups
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova

# manova for each paired groups
t1$cal_manova(manova_all = FALSE)
t1$res_manova

# PERMDISP is implemented to check multivariate homogeneity of groups dispersions (variances)
# for the whole comparison and for each paired groups
t1$cal_betadisper()
t1$res_betadisper


##########################################################################
###  Model-based class ###

# the differential test result is stored in the object$res_diff
# run metastat example
# metastat analysis at Genus level
t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Unified", taxa_level = "Genus")
# t1$res_diff is the differential test result
# t1$res_abund is the group abundance

# the metastat can run the comparisons for each paired group
# the user should use select_group to select the required pair

# select_group should be one of groups in t1$res_diff$Comparison
# FC_f vs FC_c differential abundance
t1$plot_diff_abund(use_number = 1:20, col = c("tomato2", "tomato4"), select_group = "FC_f - FC_c")
#  save plot with 600 dpi resolution
dev.print(tiff, "FC_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# PL_f vs PL_c differential abundance
t1$plot_diff_abund(use_number = 1:20, col = c("seagreen2", "seagreen4"), select_group = "PL_f - PL_c")
#  save plot with 600 dpi resolution
dev.print(tiff, "PL_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# CM_f vs CM_c differential abundance
t1$plot_diff_abund(use_number = 1:20, col = c("steelblue2", "steelblue4"), select_group = "CM_f - CM_c")
#  save plot with 600 dpi resolution
dev.print(tiff, "CM_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# LEfSe combines the non-parametric test and linear discriminant analysis
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Type", alpha = 0.01, lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 4)
#  save plot with 600 dpi resolution
dev.print(tiff, "LEfSE_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# ANOVA method and transpose
#install.packages("agricolae")
t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Unified", taxa_level = "Family", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:10, col = c1, add_sig = T, coord_flip = F)
#  save plot with 600 dpi resolution
dev.print(tiff, "ANOVA_diff_abundance.tiff", compression = "lzw", res=800, height=15, width=35, units="cm")


# Try to run those examples
# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset, method = "KW", group = "Unified", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:10, col = c1, add_sig = T, coord_flip = F)
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
t1$cal_ordination(method = "RDA", taxa_level = "Genus")
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
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

# t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "Group", plot_shape = "Type", col = c2, plot_type = c("point", "ellipse"))
#  save plot with 600 dpi resolution
dev.print(tiff, "RDA_genus.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# Mantel test can be used to check significant correlations between environmental variables and distance matrix
t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
t1$res_mantel


# perform a correlation heatmap between environmental and taxa data at Phylum level
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
t1$cal_cor(use_data = "Family", cor_method = "spearman", p_adjust_method = "fdr")
# return t1$res_cor

# plot the correlation results using plot_cor function
# filter phylum that do not have at least one ***
t1$plot_cor(color_vector = c("#A50026","white","#053061"), keep_full_name = F, keep_prefix = F)
#  save plot with 600 dpi resolution
dev.print(tiff, "phylum_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# pheatmap method
# clustering heatmap; require pheatmap package
# Let's take another color pallete: 
#install.packages('pheatmap')
t1$plot_cor(pheatmap = TRUE, color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))
#  save plot with 600 dpi resolution
dev.print(tiff, "pheatmap_phylum_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


#correlations between environmental variables and taxa for different groups
# calculate correlations for different groups using parameter by_group

# first create trans_diff object as a demonstration
t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Unified", rf_taxa_level = "Genus")
# then create trans_env object
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
# use other_taxa to select taxa you need
t1$cal_cor(by_group = "Unified", use_data = "other", cor_method = "spearman", 
           p_adjust_method = "fdr", other_taxa = t2$res_diff$Taxa[1:25])
t1$plot_cor()
#  save plot with 600 dpi resolution
dev.print(tiff, "env_corr.tiff", compression = "lzw", res=600, height=15, width=40, units="cm")


# relationship between environmental factors and alpha diversity
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
# use add_abund_table parameter to add the extra data table
t1$cal_cor(add_abund_table = dataset$alpha_diversity, cor_method = "spearman")
# Let's try to use ggplot2 with clustering plot
t1$plot_cor(cluster_ggplot = "row")
#  save plot with 600 dpi resolution
dev.print(tiff, "env_div_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")
