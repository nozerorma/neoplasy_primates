####################################
##COMPONENTS STUDY FOR TRAITS GROUPS##
###BASED ON PRESENCE/ABSCENSE, TO ###
###DETERMINE CONSISTENCY FROM DIFFERENT ##
## # SOURCE DATASETS AND KNOW THE TRAITS 
##WE CAN COMPARE BETWEEN EACH OTHER DUE TO ###
#THE SPECIES THEY CONBTAIN ###########
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
#library(ComplexHeatmap)
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

#A) READ LAST TRAITS FILE and metadata

all_traits <- read.csv("~/Documents/PhD_EvoGenomics/1st_year/GenomePhenome/data/Phenome/AllTraits25-06-20.tsv", header=T, sep = "\t")

#B) SELECT DESIRED SPECIES (sequenced_ones) WITH GROUPS IN PHYLOGENY

seq_species <- read.csv("/home/alejandro1395/Documents/PhD_EvoGenomics/1st_year/GenomePhenome/data/Phenome/names.txt",
                        sep = "\t", header = TRUE)
total_seq_species <- as.character(seq_species$Trait_names[!is.na(seq_species$Trait_names)])

sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% total_seq_species),]

#ADD GROUP INFO
groups_list <- strsplit(as.character(all_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups

all_traits <- cbind(all_traits, groups)


#C) SELECT FILE WITH CANCER TRAITS SPECIES COVERED

cancer_species <- read.csv("../data/Phenome/cancer_species.txt", header=T, sep = "\t")

#fill phylogenetic information

for (i in 1:length(cancer_species$Cancer_species)){
  if(as.character(cancer_species$Cancer_species[i]) %in% all_traits$SpeciesBROAD){
  cancer_species$Group[i] <- as.character(all_traits$groups[as.character(all_traits$SpeciesBROAD) == as.character(cancer_species$Cancer_species[i])])
  }
  else{
    cancer_species$Group[i] <- NA
  }}

View(cancer_species)


#cancer traits
cancer_species_traits <- read.csv("../data/Phenome/cancer_traits.csv", header=T, sep = "\t")

#species names
species_names <- read.csv("../data/Phenome/names.txt", header=T, sep = "\t")

#tree
tree <- read.newick("../data/Tree/primates.newick")

tree_species <- c()

###################################
#ANALYSIS WITH ONLY PUBLIC GENOMES#
###################################

overlapping_species <-cancer_species$Cancer_species[which(cancer_species$UCSC_alignments == 1)]
overlapping_families <- cancer_species$Group[which(cancer_species$UCSC_alignments == 1)]
SelectedSpecies <- cancer_species_traits[which(cancer_species_traits$Species %in% overlapping_species),]


for (i in 1:length(SelectedSpecies$Species)){
  if (SelectedSpecies$Species[i] %in% species_names$Cancer_species){
    tree_species[i] <- as.character(species_names$Tree_species[which(species_names$Cancer_species == as.character(SelectedSpecies$Species[i]))])
  }}


pruned_primates.tree<-drop.tip(tree,tree$tip.label[-match(tree_species, tree$tip.label)])


Malignancy <- SelectedSpecies$PropMalignant
names(Malignancy) <- SelectedSpecies$Species

obj<-contMap(pruned_primates.tree,Malignancy,plot=FALSE)
obj<-setMap(obj,c("#53EC2D", "#F1F97C", 
                  "white", "#F08431", "#F03131"))
plot(obj,fsize=c(1),lwd=c(5),leg.txt="Malignancy proportion")
lower_bound <- quantile(Malignancy, 0.25)
upper_bound <- quantile(Malignancy, 0.75)
outlier_ind <- which(Malignancy < lower_bound | Malignancy > upper_bound)


Benign <- SelectedSpecies$PropBenign
names(Benign) <- SelectedSpecies$Species

obj<-contMap(pruned_primates.tree,Benign,plot=FALSE)
obj<-setMap(obj,rev(c("#53EC2D", "#F1F97C", 
                  "white", "#F08431", "#F03131")))
plot(obj,fsize=c(1),lwd=c(5),leg.txt="Benign proportion")

lower_bound <- quantile(Benign, 0.25)
upper_bound <- quantile(Benign, 0.75)
outlier_ind <- which(Benign < lower_bound | Benign > upper_bound)




Neoplasy<- SelectedSpecies$NeoplasiaRate
names(Neoplasy) <- SelectedSpecies$Species

obj<-contMap(pruned_primates.tree,Neoplasy,plot=FALSE)
obj<-setMap(obj,c("#53EC2D", "#F1F97C", 
                  "white", "#F08431", "#F03131"))
plot(obj,fsize=c(1),lwd=c(5),leg.txt="Neoplasy rate")
lower_bound <- quantile(Neoplasy, 0.25)
upper_bound <- quantile(Neoplasy, 0.75)
outlier_ind <- which(Neoplasy > upper_bound)
pf <- twoSampleFactor(as.matrix(Benign), pruned_primates.tree, nfactors=1)
pp <- pf.taxa(pf,taxonomy,factor=1)
summary(pf,taxonomy,factor=1)



#BAD SPECIES --> Macaca fascicularis, Macaca mulatta, Porpithecus y microcebus murinus
#GOOD species --> Colobus angolensis, sabaeus, Pan paniscus y troglodytes, Cebus capucinus,
#sAIMIRI BOLIVIENSIS)

#let's try phylofactor'
library(phylofactor)
library(ggtree)
taxonomy <- as.data.frame(cbind(as.character(overlapping_species),
                          as.character(overlapping_families)))
colnames(taxonomy) <- c("species", "family")
pf <- twoSampleFactor(as.matrix(Malignancy), pruned_primates.tree, nfactors=2)
pp <- pf.taxa(pf,taxonomy,factor=2)
  
#No me gusta, te separa por clados NOOOOOOOOO

data <- kmeans(Benign, 3)
clusters <- data$cluster
  
library(ggtree)
p <- ggtree(pruned_primates.tree)
data_clusters <- as.data.frame(cbind(names(clusters),
                    clusters))
colnames(data_clusters) <- c("species", "cluster")
p +  geom_tippoint(aes(color=data_clusters$cluster))












#Transform data for plot

rownames(cancer_species) <- cancer_species$Cancer_species

cancer_species_df <- as.data.frame(melt(cancer_species))

View(cancer_species)
cancer_species_counts <- aggregate(as.numeric(as.character(cancer_species_df$value)),
          by=list(cancer_species_df$Group, 
                  cancer_species_df$variable),
          FUN=sum)
colnames(cancer_species_counts) <- c("Family", "Source", "Count")

primates_palette <- readRDS("~/Documents/PhD_EvoGenomics/1st_year/GenomePhenome/src/primate_families_palette.rds")

#How is in the case of the species we have

table_total_counts <- table(cancer_species$Group)
total_df <- as.data.frame(cbind(names(table_total_counts), 
                    rep("total", length(table_total_counts)),
                    table_total_counts), row.names = FALSE)
colnames(total_df) <- c("Family", "Source", "Count")

total_coverage_df <- rbind(cancer_species_counts,
                           total_df)
#How they cover our dataset
ggplot(total_coverage_df) +
  geom_bar(aes(x=Source,
               fill=Family,
               y=as.numeric(as.character(Count))),
           stat="identity") + scale_fill_manual(values = primates_palette) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 25,
                                                     size=10)) +
  ylab("Count species")




###############################################
##Now let's see what happens with the traits###
###############################################

cancer_species_traits <- read.csv("../data/Phenome/cancer_traits.csv", header=T, sep = "\t")


#total variance analysis

for (i in 1:length(cancer_species_traits$Species)){
  if(as.character(cancer_species_traits$Species[i]) %in% all_traits$SpeciesBROAD){
    cancer_species_traits$Group[i] <- as.character(all_traits$groups[as.character(all_traits$SpeciesBROAD) == as.character(cancer_species_traits$Species[i])])
  }
  else{
    cancer_species_traits$Group[i] <- NA
  }}

View(cancer_species_traits)

#all are numeric variables

GrouppedCancer_traits_alldf <- GrouppedCancer_traits[,3:4]
ContinousTraitsMeans <- GrouppedCancer_traits_alldf %>%
  group_by(variable) %>%
  summarize_all(mean, na.rm = TRUE)

ContinousTraitsSds<- GrouppedCancer_traits_alldf %>%
  group_by(variable) %>%
  summarize_all(sd, na.rm = TRUE)


AllStatsContinous <- merge(ContinousTraitsMeans, ContinousTraitsSds, by=c("variable"))
colnames(AllStatsContinous) <- c("Trait", "Mean", "Sd")
AllStatsContinous$coef_var <- AllStatsContinous$Sd/AllStatsContinous$Mean

ggplot(data=AllStatsContinous) + 
  geom_bar(aes(x=Trait, 
               y=as.numeric(as.character(coef_var))), 
           stat="identity",
           colour="black") +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.5)) +
  ylab("Coeficient of variance")




#by each family information
GrouppedCancer_traits <- melt(cancer_species_traits)
GrouppedCancer_traits_df <- GrouppedCancer_traits[,2:4]
ContinousTraitsMeans_clades <- GrouppedCancer_traits_df %>%
  group_by(Group, variable) %>%
  summarize_all(mean, na.rm = TRUE)

ContinousTraitsSds_clades <- GrouppedCancer_traits_df %>%
  group_by(Group, variable) %>%
  summarize_all(sd, na.rm = TRUE)

GrouppedCancer_traits_df[which(GrouppedCancer_traits_df$Group == "PLA_Aotidae"),]

AllStatsContinous_clades <- merge(ContinousTraitsMeans_clades, ContinousTraitsSds_clades, by=c("Group", "variable"))
colnames(AllStatsContinous_clades) <- c("Family", "Trait", "Mean", "Sd")
AllStatsContinous_clades$coef_var <- AllStatsContinous_clades$Sd/AllStatsContinous_clades$Mean*100

ggplot(data=AllStatsContinous_clades) + 
  geom_bar(aes(x=Family, 
               y=as.numeric(as.character(coef_var)), 
               fill=Family), 
           stat="identity",
           colour="black") +
  scale_fill_manual(values = primates_palette) + 
  theme(legend.title = element_blank()) +
  facet_wrap(~Trait) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Coeficient of variance")


###################################
#ANALYSIS WITH ONLY PUBLIC GENOMES#
###################################

overlapping_species <-cancer_species$Cancer_species[which(cancer_species$UCSC_alignments == 1)]
SelectedSpecies <- GrouppedCancer_traits[which(GrouppedCancer_traits$Species %in% overlapping_species),]


#all species

GrouppedCancer_traits_alldf <- SelectedSpecies[,3:4]
ContinousTraitsMeans <- GrouppedCancer_traits_alldf %>%
  group_by(variable) %>%
  summarize_all(mean, na.rm = TRUE)

ContinousTraitsSds<- GrouppedCancer_traits_alldf %>%
  group_by(variable) %>%
  summarize_all(sd, na.rm = TRUE)


AllStatsContinous <- merge(ContinousTraitsMeans, ContinousTraitsSds, by=c("variable"))
colnames(AllStatsContinous) <- c("Trait", "Mean", "Sd")
AllStatsContinous$coef_var <- AllStatsContinous$Sd/AllStatsContinous$Mean

ggplot(data=AllStatsContinous) + 
  geom_bar(aes(x=Trait, 
               y=as.numeric(as.character(coef_var))), 
           stat="identity",
           colour="black") +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.5)) +
  ylab("Coeficient of variance")


#by clade
GrouppedCancer_traits_df <- SelectedSpecies[,2:4]
ContinousTraitsMeans_clades <- GrouppedCancer_traits_df %>%
  group_by(Group, variable) %>%
  summarize_all(mean, na.rm = TRUE)

ContinousTraitsSds_clades <- GrouppedCancer_traits_df %>%
  group_by(Group, variable) %>%
  summarize_all(sd, na.rm = TRUE)

GrouppedCancer_traits_df[which(GrouppedCancer_traits_df$Group == "PLA_Aotidae"),]

AllStatsContinous_clades <- merge(ContinousTraitsMeans_clades, ContinousTraitsSds_clades, by=c("Group", "variable"))
colnames(AllStatsContinous_clades) <- c("Family", "Trait", "Mean", "Sd")
AllStatsContinous_clades$coef_var <- AllStatsContinous_clades$Sd/AllStatsContinous_clades$Mean

ggplot(data=AllStatsContinous_clades) + 
  geom_bar(aes(x=Family, 
               y=as.numeric(as.character(coef_var)), 
               fill=Family), 
           stat="identity",
           colour="black") +
  scale_fill_manual(values = primates_palette) + 
  theme(legend.title = element_blank()) +
  facet_wrap(~Trait) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Coeficient of variance")


overlapping_traits <- cancer_species_traits[which(cancer_species_traits$Species %in% overlapping_species),]


var_quanti <- PCA(overlapping_traits[,2:16], 
                  scale.unit = TRUE)
fviz_pca_var(var_quanti, col.var="contrib",
             repel = TRUE)
var_quanti_cor <- as.data.frame(var_quanti$var$cor)
library(Hmisc)
library(corrplot)
res <- cor(overlapping_traits[,2:16],
           use="pairwise.complete.obs")
corrplot(res,
         tl.col = "black")

#First axis mainly driven by correlations with the Total nº of records
