# First approach to non-taxonomic phenotipic distribution
## The idea is to group the species given their branch lengths instead of the 
## archeotypical taxon view, which is somewhat biassed.

set.seed(22)
library(ape)
library(RRphylo)

# Tree import

tree <- read.newick("Data/233-GENOMES/science.abn7829_data_s4.nex.tree")


# Read csv for neoplasia trait

cancer_species_traits <- read.csv("Data/Neoplasia/species360_primates_neoplasia_20230519.pruned.csv", header=T, sep = ",")

# Normalize name using binomial syntax
cancer_species_traits$SPECIES_BINOMIAL <- sub(" ", "_", cancer_species_traits$species)
cancer_species_traits$SPECIES_BINOMIAL -> speciesbin
speciesbin -> names(speciesbin)

pruned_tree <- treedataMatch(tree=tree,y=speciesbin)
removed_tree <- pruned_tree$removed.from.tree
removed_phenotype <- pruned_tree$removed.from.y

# Are the missing ones really missing?

# Extract the prefix (before the underscore) for each element in removed_phenotype
prefixes <- sapply(strsplit(removed_phenotype, "_"), `[`, 1)
suffixes <- sapply(strsplit(removed_phenotype, "_"), `[`,2)

# Convert both vectors to lowercase for case-insensitive matching
prefixes_lower <- tolower(prefixes)
suffixes_lower <- tolower(suffixes)
tree_tips_lower <- tolower(tree$tip.label)

# Identify any tree tip label that starts with the extracted prefix
matching_tips <- sapply(prefixes_lower, function(prefix) {
  grep(paste0("^", prefix), tree_tips_lower)
})
matching_tips_suf <- sapply(suffixes_lower, function(suffix) {
  grep(paste0(suffix), tree_tips_lower)
})

# Flatten the list to a vector
matched_indices <- unique(unlist(matching_tips))


# Extract the labels that matched
matched_tree_tips <- tree$tip.label[matched_indices]
matched_tree_tips_suf <- tree$tip.label[matched_indices_suf]

# export pruned tree
write.tree(pruned_tree$tree, file="Data/233-GENOMES/science.abn7829_data_s4.nex.tree.pruned")


# Okay now let's work with the tree
## Visualization

ggtree(pruned_tree$tree) + 
  geom_tiplab() + 
  geom_treescale(x=0, fontsize=3) +
  theme(legend.position="bottom")

# Generate tree cuts at given ages

at30 <- cutPhylo(pruned_tree$tree,age=30)
at20 <- cutPhylo(pruned_tree$tree,age=20)


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
