### data analysis with microeco package - agricultural waste paper - Guilherme Martins (16/07/2022)

# load the package
library(microeco)

# Get and set work directory
getwd()
setwd("C:/Users/glmartins/OneDrive/Doutorado/Colaborações/Experimento glifosato")
path <- "C:/Users/glmartins/OneDrive/Doutorado/Colaborações/Experimento glifosato"
list.files(path)

# load the example data; 16S rRNA gene amplicon sequencing dataset
# metadata table
sample_info_ITS <- read.table("sample_info_ITS.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# feature table
otu_table_ITS <- read.table("otu_table_ITS.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# taxonomic assignment table
taxonomy_table_ITS <- read.table("taxonomy_table_ITS.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# load the environmental data table if it is not in sample table
env_data_ITS <- read.table("env_data_ITS.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# use pipe operator in magrittr package
library(magrittr) # needs to be run every time you start R and want to use %>%

# set.seed is used to fix the random number generation to make the results repeatable
set.seed(123)

# make the plotting background same with the tutorial
library(ggplot2)
theme_set(theme_bw())

# Make sure that the data types of sample_table, otu_table and tax_table are all data.frame format as the following part shows
class(otu_table_ITS)
otu_table_ITS[1:5, 1:5]

class(taxonomy_table_ITS)
taxonomy_table_ITS[1:5, 1:3]


#######################################################################
            ### Organize data for further analysis ###


# Remove NAs to proceed with diversity analysis
# make the taxonomic information unified, very important
taxonomy_table_ITS %<>% tidy_taxonomy

# make sure that the rownames of sample information table are sample names
class(sample_info_ITS)
sample_info_ITS[1:5, ]

# make sure that the environmental data are stored as env data
class(env_data_ITS)
env_data_ITS[1:5, 1:5]

# In R6 class, '$new' is the original method used to create a new object of class
# Let's create a microtable object with more information
dataset <- microtable$new(sample_table = sample_info_ITS, otu_table = otu_table_ITS, tax_table = taxonomy_table_ITS)
dataset

# make the OTU and sample information consistent across all files in the dataset object
dataset$tidy_dataset()
print(dataset)

# let's use sample_sums() to check the sequence numbers in each sample
dataset$sample_sums() %>% range

# Rarefy data to reduce bias in diversity index
dataset$rarefy_samples(sample.size = 1359)

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

#######################################################################
### Calculate alpha and beta diversity ###

### Alpha diversity ###

# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = FALSE)

# return dataset$alpha_diversity
class(dataset$alpha_diversity)

# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")

# calculate beta diversity

# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")


################################################################################
### Composition-based class ###

# taxa_abund list in the object of microtable class must be first calculated
# create trans_abund object


# reorder groups
dataset$sample_table$Soil %<>% factor(., levels = c("NF", "CT", "CTN", "NT", "NTN"))
dataset$sample_table$Dilution %<>% factor(., levels = c("CS", "D1", "D3"))
dataset$sample_table$Code %<>% factor(., levels = c("NF", "NF_D1", "NF_D3", "CT", "CT_D1", "CT_D3", "CTN", "CTN_D1", "CTN_D3", "NT", "NT_D1", "NT_D3", "NTN", "NTN_D1", "NTN_D3"))
str(dataset$sample_table)


##########################################################################
### Define color palette ###
c1 <- c("seagreen4", "seagreen3", "seagreen2",
        "mediumpurple4", "mediumpurple3", "mediumpurple2",
        "orange4", "orange3", "orange2",
        "steelblue4", "steelblue3", "steelblue2",
        "tomato4", "tomato3", "tomato2")

c2 <- c("seagreen3", "mediumpurple3", "orange3","steelblue3", "tomato3")

##########################################################################


# use 10 Phyla with the highest abundance in the dataset.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 10)
# Relative abundance with two facets at Phylum level
# require package ggh4x, first run install.packages("ggh4x") if not installed
#install.packages("ggh4x")
t1$plot_bar(others_color = "grey70", facet = "Soil", facet2 = "Dilution", xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)
# save plot with 600 dpi resolution
dev.print(tiff, "relative_abundance_facet2.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# show the heatmap with the high abundant genera
# show 30 taxa at Genus level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 30)
t1$plot_heatmap(facet = "Soil", xtext_keep = FALSE, withmargin = TRUE)
#  save plot with 600 dpi resolution
dev.print(tiff, "heatmap_abundance4.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# The trans_venn class is used for venn analysis, i.e. shared and unique taxa
# To analyze the unique and shared OTUs of groups, we first merge samples according to the "Group" column of sample_table

# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Dilution")
# dataset1 is a new microtable object
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn() # The integer is OTU number # The percentage data is the sequence number/total sequence number
#  save plot with 600 dpi resolution
dev.print(tiff, "venn_diagram.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# transform the results of venn plot to the traditional feature-sample table, that is, another object of microtable class
dataset1 <- dataset$merge_samples(use_group = "Dilution")
t1 <- trans_venn$new(dataset1)

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)


# Venn diagram by group
dataset1 <- dataset$merge_samples(use_group = "Code")
# dataset1 is a new microtable object
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn(petal_plot = TRUE) # The integer is OTU number # The percentage data is the sequence number/total sequence number
#  save plot with 600 dpi resolution
dev.print(tiff, "venn_diagram3.tiff", compression = "lzw", res=600, height=15, width=15, units="cm")


################################################################################
### Diversity-based class ###

# The data_alpha is used for the following differential test and plotting
t1 <- trans_alpha$new(dataset = dataset, group = "Code")
# return t1$data_stat
t1$data_stat[1:5, ]

# test the differences among groups using:
# Kruskal-Wallis Rank Sum Test (overall test when groups > 2)
t1$cal_diff(method = "KW")
# return t1$res_diff
t1$res_diff[1:5, ]

# Dunn's Kruskal-Wallis Multiple Comparisons (for paired groups when groups > 2)
#install.packages("FSA")
#install.packages("agricolae")

t1$cal_diff(method = "KW_dunn")
# return t1$res_diff
t1$res_diff[1:5, ]

# ANOVA with multiple comparisons
t1$cal_diff(method = "anova")
# return t1$res_diff
t1$res_diff[1:5, ]

# plot the mean and se of alpha diversity for each group, and add the anova result
t1$cal_diff(method = "anova")
t1$plot_alpha(measure = "Chao1", color_values = c1)
t1$plot_alpha(measure = "ACE", color_values = c1)
t1$plot_alpha(measure = "Simpson", color_values = c1)
t1$plot_alpha(measure = "Shannon", color_values = c1)
t1$plot_alpha(measure = "Observed", color_values = c1)


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
t1$plot_clustering(group = "Code", replace_name = c("Soil", "Dilution"))
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
t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Code", taxa_level = "Genus", use_number = 6)
# t1$res_diff is the differential test result
# t1$res_abund is the group abundance

# the metastat can run the comparisons for each paired group
# the user should use select_group to select the required pair

# select_group should be one of groups in t1$res_diff$Comparison
# NF vs NF_D3 differential abundance
fun_FN <- F_NF <- t1$plot_diff_abund(
  use_number = 1:10,
  col = c("seagreen2", "seagreen4"),
  add_sig = T,
  coord_flip = T,
  select_group = "NF - NF_D3",
  text_y_size = 10,
  keep_prefix = F,
  simplify_names = T,
  barwidth = 0.6
  )

fun_FN <- F_NF + rremove("xlab")
fun_FN

#  save plot with 600 dpi resolution
dev.print(tiff, "NF_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# CF vs CF_D3 differential abundance
fun_CF <- F_CF <- t1$plot_diff_abund(use_number = 1:10,
                           col = c("mediumpurple2", "mediumpurple4"),
                           add_sig = T,
                           coord_flip = T,
                           select_group = "CF - CF_D3",
                           text_y_size = 10,
                           keep_prefix = F,
                           simplify_names = TRUE,
                           barwidth = 0.6
                           )

fun_CF <- F_CF + rremove("xlab")
fun_CF

#  save plot with 600 dpi resolution
dev.print(tiff, "CF_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# CFN vs CFN_D3 differential abundance
t1$plot_diff_abund(use_number = 1:10,
                   col = c("orange2", "orange4"),
                   add_sig = T,
                   coord_flip = F,
                   select_group = "CFN - CFN_D3",
                   text_y_size = 10,
                   keep_prefix = F,
                   keep_prefix = F,
                   simplify_names = TRUE,
                   barwidth = 0.6
                   )

#  save plot with 600 dpi resolution
dev.print(tiff, "CFN_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")

# NT vs NT_D3 differential abundance
fun_NT <- F_NT <- t1$plot_diff_abund(use_number = 1:10,
                           col = c("steelblue2", "steelblue4"),
                           add_sig = T,
                           coord_flip = T,
                           select_group = "NT - NT_D3",
                           text_y_size = 10,
                           keep_prefix = F,
                           simplify_names = TRUE,
                           barwidth = 0.6
                           )

fun_NF <- F_NT + rremove("xlab")
fun_NF

#  save plot with 600 dpi resolution
dev.print(tiff, "NT_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")

# NTN vs NTN_D3 differential abundance
t1$plot_diff_abund(use_number = 1:10,
                   col = c( "tomato2", "tomato4"),
                   add_sig = T,
                   coord_flip = F,
                   select_group = "NTN - NTN_D3"
                   )

#  save plot with 600 dpi resolution
dev.print(tiff, "NTN_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")

# Combine all diff abundance analysis in one figure
library("ggpubr")
fung_diff <- ggarrange(fun_FN, fun_NT, fun_CF, 
          labels = c("NF", "NT", "CF"),
          common.legend = FALSE,
          legend = "none",
          align = "hv",
          hjust = -10,
          vjust = 1.3,
          widths = c(1,1,1),
          heights = c(1,1,1),
          ncol = 3, nrow = 1)

fung_diff

#  save plot with 800 dpi resolution
dev.print(tiff, "Fun_diff_abundance.tiff", compression = "lzw", res=800, height=15, width=40, units="cm")


# LEfSe combines the non-parametric test and linear discriminant analysis
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Dilution", alpha = 0.01, lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 4)
#  save plot with 600 dpi resolution
dev.print(tiff, "LEfSE_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# ANOVA method and transpose
#install.packages("agricolae")
t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Code", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:10, col = c1, add_sig = T, coord_flip = F)
#  save plot with 600 dpi resolution
dev.print(tiff, "ANOVA_diff_abundance.tiff", compression = "lzw", res=600, height=20, width=50, units="cm")


# Try to run those examples
# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset, method = "KW", group = "Code", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:10, col = c1, add_sig = T, coord_flip = T)
#  save plot with 600 dpi resolution
dev.print(tiff, "KW_diff_abundance.tiff", compression = "lzw", res=600, height=25, width=20, units="cm")


#################################################################################
### Explainable class ###

# If there may be some NA in the user's env data,add complete_na = TRUE when creating the trans_env object


# add_data is used to add the environmental data
t1 <- trans_env$new(dataset = dataset, add_data = env_data_ITS)

# show the autocorrelations among variables
# use group parameter to show the distributions of variables and the autocorrelations across groups
#install.packages("GGally")
t1$cal_autocor(group = "Soil", col = c2)
#  save plot with 600 dpi resolution
dev.print(tiff, "autocorr2.tiff", compression = "lzw", res=600, height=40, width=40, units="cm")


# RDA at Genus level
t1$cal_ordination(method = "RDA", taxa_level = "Genus")
t1$res_ordination
t1$res_ordination_R2
t1$res_ordination_terms
t1$res_ordination_axis

#r.squared adj.r.squared 
#0.6277424     0.5573153

#Df Variance       F Pr(>F)    
#B_glucosidase     1    29362  9.8563  0.001 ***
#Acid_phosphatase  1    50342 16.8990  0.001 ***
#pH                1    29104  9.7697  0.001 ***
#M.O.              1    38282 12.8504  0.001 ***
#P                 1    12289  4.1251  0.004 ** 
#K                 1     4136  1.3883  0.205    
#C.CO2             1    22356  7.5046  0.001 ***
#Residual         37   110223      

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

# t1$res_rda_trans is the transformed result for plotting
RDA_ITS <- t1$plot_ordination(plot_color = "Soil", plot_shape = "Dilution", col = c2, plot_type = c("point", "ellipse"), ellipse_chull_fill = F)
g2 <- RDA_ITS + scale_shape_manual(values=c(15, 16, 17, 18)) + geom_text_repel()
g20 <- g2 + ggtitle("Fungal community")

g20

#  save plot with 600 dpi resolution
dev.print(tiff, "RDA_genus_ITS.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# Mantel test can be used to check significant correlations between environmental variables and distance matrix
t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
t1$res_mantel


# perform a correlation heatmap between environmental and taxa data at Phylum level
t1 <- trans_env$new(dataset = dataset, add_data = env_data_ITS)
t1$cal_cor(use_data = "Genus", cor_method = "spearman", p_adjust_method = "fdr")
# return t1$res_cor

# plot the correlation results using plot_cor function
# filter phylum that do not have at least one ***
t1$plot_cor()
#  save plot with 600 dpi resolution
dev.print(tiff, "phylum_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# pheatmap method
# clustering heatmap; require pheatmap package
# Let's take another color pallete: 
install.packages('pheatmap')
t1$plot_cor(pheatmap = TRUE, color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))
#  save plot with 600 dpi resolution
dev.print(tiff, "pheatmap_phylum_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


#correlations between environmental variables and taxa for different groups
# calculate correlations for different groups using parameter by_group

# first create trans_diff object as a demonstration
t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Code", rf_taxa_level = "Genus")
# then create trans_env object
t1 <- trans_env$new(dataset = dataset, add_data = env_data_ITS)
# use other_taxa to select taxa you need
t1$cal_cor(by_group = "Code", use_data = "other", cor_method = "spearman", 
           p_adjust_method = "fdr", other_taxa = t2$res_diff$Taxa[1:25])
t1$plot_cor()
#  save plot with 600 dpi resolution
dev.print(tiff, "env_corr.tiff", compression = "lzw", res=600, height=15, width=40, units="cm")


# relationship between environmental factors and alpha diversity
t1 <- trans_env$new(dataset = dataset, add_data = env_data_ITS)
# use add_abund_table parameter to add the extra data table
t1$cal_cor(add_abund_table = dataset$alpha_diversity, cor_method = "spearman")
# Let's try to use ggplot2 with clustering plot
#install.packages("factoextra")
#install.packages("aplot")
t1$plot_cor(cluster_ggplot = "row")
#  save plot with 600 dpi resolution
dev.print(tiff, "env_div_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# FUNGuild and Fungal Traits function
install.packages('igraph')

# create microtable object
meco_fungi <- microtable$new(sample_table = sample_info_ITS, otu_table = otu_table_ITS, tax_table = taxonomy_table_ITS)
meco_fungi

# remove the taxa not assigned in the Kingdom "k__Fungi"
meco_fungi$tax_table %<>% base::subset(Kingdom == "k__Fungi")

# use tidy_dataset() to make OTUs and samples information consistent across files
meco_fungi$tidy_dataset()

# create trans_network object
t1 <- trans_network$new(dataset = meco_fungi, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0.000001, cor_method = "spearman", node_label = "name")
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# add modules
#t1$cal_module()

# convert module info to microtable object
#meco_module <- t1$trans_comm(use_col = "module")

# create trans_func object
t2 <- trans_func$new(meco_fungi)

# identify species traits, automatically select database for prokaryotes or fungi
# fungi_database = "FungalTraits" for the FungalTraits database
t2$cal_spe_func(fungi_database = "FUNGuild")

# calculate abundance-unweighted functional redundancy of each trait for each network module
t2$cal_spe_func_perc(abundance_weighted = FALSE)
t2$res_spe_func_perc

# Extract FUN GUILD table
funguild <- t2$res_spe_func_perc
funguild2 <- funguild[, -c(4, 18, 20, 21, 22)]

# Export table
write.table(funguild2, file = "funguild_table.txt", sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)


# plot the functional redundancy of network modules
t2$plot_spe_func_perc()

