# First approach to non-taxonomic phenotipic distribution
## The idea is to group the species given their branch lengths instead of the 
## archeotypical taxon view, which is somewhat biassed.

# Tree import

# NOT RUN {
library(ape)
library(RRphylo)

set.seed(22)
rtree(100)->tree
30->age

at30 <- cutPhylo(tree,age=30)
at20 <- cutPhylo(tree,age=20)
cutPhylo(tree,node=30)->t2
cutPhylo(tree,node=30,keep.lineage=TRUE)->t2
# }

library(ape)

ggtree(tree) + geom_tippoint(aes(color=factor(at30)))

df <- data.frame(label = tree$tip.label, cluster = clusters_wardd2)
ggtree(tree) + geom_tippoint(aes(color=factor(cluster)), data=df) + scale_color_discrete(name="Cluster")

tip_labels <- tree$tip.label
reordered_clusters <- clusters_wardd2[tip_labels]

ggtree(tree) + geom_tippoint(aes(color=factor(reordered_clusters))) + scale_color_discrete(name="Cluster")

######################################
####All species trait distribution####
###Based solely on neoplasia ratio####
######################################







library(ggplot2)
library(gridExtra)
library(tibble)
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(dendextend)
library(ComplexHeatmap)
library("factoextra")
library("circlize")
library("RColorBrewer")
library(ggnewscale)
library(scales)
library(cowplot)
library(ggdendro)
library(viridis)
library(cowplot)
library(FactoMineR)
library(missMDA)
library(psych)
library(network)
library(dendextend)
library(ape)
library(randomcoloR)
library(phytools)

# Read csv for neoplasia trait

cancer_species_traits <- read.csv("/home/miguel/TFM/Master_projects/NEOPLASY_PRIMATES/Data/Neoplasia/species360_primates_neoplasia_20230519.csv", header=T, sep = "\t")

# Normalize name using binomial syntax
cancer_species_traits$SPECIES_BINOMIAL <- sub(" ", "_", cancer_species_traits$species)

# Tree input
tree <- read.newick("./Data/233-GENOMES/science.abn7829_data_s3.nw.tree")

# Histogramize (?)
hist(cancer_species_traits$neoplasia_prevalence)

# What's done here?
library(e1071)
skewness(cancer_species_traits$neoplasia_prevalence)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

tree_species <- c()
tree_species_names <- c()
tree_species_names$Tree_species<- tree[["tip.label"]]
view(tree_species_names)

###################################
#ANALYSIS WITH ONLY TREE SPECIES#
###################################

for (i in 1:length(cancer_species_traits$SPECIES_BINOMIAL)){
    tree_species[i] <- as.character(tree_species_names$Tree_species[which(tree_species_names$SPECIES_BINOMIAL == as.character(cancer_species_traits$SPECIES_BINOMIAL[i]))])
  }

pruned_primates.tree<-drop.tip(tree,tree$tip.label[-match(tree_species[!is.na(tree_species)], 
                                                          tree$tip.label)])

for (i in 1:length(SelectedSpecies$Species)){
    tree_species[i] <- as.character(tree_species_names$Tree_species[which(species_names$Cancer_species == as.character(SelectedSpecies$Species[i]))])
  }


pruned_primates.tree<-drop.tip(tree,tree$tip.label[-match(cancer_species_traits$SPECIES_BINOMIAL, tree$tip.label)])
