##################
###CANCER TRAITS##
##################

############################

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

cancer_species_traits <- read.csv("../data/Phenome/cancer_traits.csv", header=T, sep = "\t")
for (i in 1:length(cancer_species_traits$Species)){
  if(as.character(cancer_species_traits$Species[i]) %in% all_traits$SpeciesBROAD){
    cancer_species_traits$Group[i] <- as.character(all_traits$groups[as.character(all_traits$SpeciesBROAD) == as.character(cancer_species_traits$Species[i])])
  }
  else{
    cancer_species_traits$Group[i] <- NA
  }}
######
#TREE

#fill phylogenetic information

primates_palette <- readRDS("~/Documents/PhD_EvoGenomics/1st_year/GenomePhenome/src/primate_families_palette.rds")
#TREE WITH INFORMATION

#tree
tree <- read.newick("../data/Tree/primates.newick")
#species names
species_names <- read.csv("../data/Phenome/names.txt", header=T, sep = "\t")

cancer_filtered_traits <- cancer_species_traits %>% 
  group_by(Species, Group) %>%
  summarise_all(funs(mean(., na.rm=TRUE)))

melted_dataframe_traits <- melt(cancer_filtered_traits, id.vars = c("Group", "Species"))
colnames(melted_dataframe_traits) <- c("Family", "Species", "Trait", "value")


is_outlier <- function(x) {
  return(x < quantile(x, 0.25,na.rm=TRUE) - 1.5 * IQR(x, na.rm = TRUE) | x > quantile(x, 0.75,na.rm=TRUE) + 1.5 * IQR(x, na.rm=TRUE))
}

#all together
melted_dataframe_traits <- melted_dataframe_traits %>%
  dplyr::group_by(Trait) %>% 
  mutate(label = case_when(
    value < quantile(value, 0.25, na.rm=TRUE) ~ "low_extreme", 
    value > quantile(value, 0.75, na.rm=TRUE) ~ "high_extreme",
    TRUE ~ "normal")) %>%
  mutate(outlier = ifelse(is_outlier(value), 
                          as.character(Species), 
                          NA))



#x family
melted_dataframe_traits <- melted_dataframe_traits %>%
  dplyr::group_by(Family, Trait) %>% 
  mutate(label_fam = case_when(
    value < quantile(value, 0.25, na.rm=TRUE) ~ "low_extreme", 
    value > quantile(value, 0.75, na.rm=TRUE) ~ "high_extreme",
    TRUE ~ "normal")
  )


tree_species <- c()
for (i in 1:length(cancer_filtered_traits$Species)){
  if (cancer_filtered_traits$Species[i] %in% species_names$Cancer_species){
    tree_species[i] <- as.character(species_names$Tree_species[which(species_names$Cancer_species == as.character(cancer_filtered_traits$Species[i]))])
  }}

pruned_primates.tree<-drop.tip(tree,tree$tip.label[-match(tree_species[!is.na(tree_species)], 
                                                          tree$tip.label)])



#######################################
#Species available from PUBLIC GENOMES#
#######################################


cancer_species <- read.csv("../data/Phenome/cancer_species.txt", header=T, sep = "\t")
overlapping_species_public <-cancer_species$Cancer_species[which(cancer_species$UCSC_alignments == 1)]
overlapping_species_200primates <-cancer_species$Cancer_species[which(cancer_species$X200PrimatesProject == 1)]
public_dataframe_species <- melted_dataframe_traits[which(melted_dataframe_traits$Species %in% overlapping_species_public),]


pruned_primates.tree<-drop.tip(tree,tree$tip.label[-match(overlapping_species_public, 
                                                          tree$tip.label)])


Malignancy_rate <- public_dataframe_species$value[public_dataframe_species$Trait == "NeoplasiaRate"]
names(Malignancy_rate) <- public_dataframe_species$Species[public_dataframe_species$Trait == "NeoplasiaRate"]

obj<-contMap(pruned_primates.tree,Malignancy_rate,plot=FALSE)
obj<-setMap(obj,c("#53EC2D", "#F1F97C", 
                  "white", "#F08431", "#F03131"))
plot(obj,fsize=c(1),lwd=c(5),leg.txt="Neoplasy rate")



###########
public_dataframe_traits <- public_dataframe_species %>%
  dplyr::group_by(Trait) %>% 
  mutate(label = case_when(
    value < quantile(value, 0.25, na.rm=TRUE) ~ "low_extreme", 
    value > quantile(value, 0.75, na.rm=TRUE) ~ "high_extreme",
    TRUE ~ "normal")) %>%
  mutate(outlier_public = ifelse(is_outlier(value), 
                                 as.character(Species), 
                                 NA))

public_dataframe_traits <- public_dataframe_traits %>%
  dplyr::group_by(Family, Trait) %>% 
  mutate(label_fam = case_when(
    value < quantile(value, 0.25, na.rm=TRUE) ~ "low_extreme", 
    value > quantile(value, 0.75, na.rm=TRUE) ~ "high_extreme",
    TRUE ~ "normal")
  )

ggplot(data=melted_dataframe_traits[which(melted_dataframe_traits$Trait == "NeoplasiaRate"),], 
       aes(x=Family, 
           y=value,
           color=Family,
           fill=Family,
           fontface="bold")) + 
  scale_color_manual(values= primates_palette) +
  scale_fill_manual(values= primates_palette) +
  geom_boxplot(outlier.shape =  "triangle",
               colour = "black",
               alpha = 0.6) +
  geom_count(colour = "black") +
  facet_wrap(~Trait, scales = "free") +
  #geom_text( fontface = "bold",
    #         position = position_jitter(height = 0.08,
     #                                   width = 0.08),
      #       size=4,
       #      aes(group = Family, 
        #         label = ifelse((label == "high_extreme") |
         #                         (label == "low_extreme") ,
          #                      as.character(Species),'')))  +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white")) +
  ylab("Trait value") +
  guides(fill = guide_legend(override.aes = list(size = 0.5)))


ggplot(data=melted_dataframe_traits[which(melted_dataframe_traits$Trait %in% c("MalignancyRate",
                                                                              "NeoplasiaRate",
                                                                             "BenignRate")),], 
       aes(x=Trait, 
           y=value,
           fontface="bold")) +
  geom_boxplot(aes(fill="white"), alpha = 0) +
  facet_wrap(~Trait, scales = "free") +
  scale_color_manual(values= primates_palette) +
  scale_fill_manual(values= "white") +
  geom_point(aes(x=Trait, 
                 y=value,
                 color = ifelse(outlier != "NA",
                                Family,
                                "white"),
                 shape = ifelse(Species %in% overlapping_species_200primates,
                               "1",
                               "2")),
             position = position_jitter(height = 0.1,
                                        width = 0.1),
             size = 2.5,
             stroke=1) +
  scale_shape_manual(values = c(8, 16)) +
  theme_minimal() + theme(axis.title.x  = element_blank(),
                          axis.text.x = element_blank(),
                          legend.position = "none")

#TRAIT ALL SPECIES DISTRIBUTION


#boxplot families
ggplot(data=melted_dataframe_traits[which(melted_dataframe_traits$Trait %in% c("MalignancyRate",
                                                                               "NeoplasiaRate",
                                                                               "BenignRate")),], 
       aes(x=Trait, 
           y=value,
           color=Family,
           fill=Family,
           fontface="bold")) + geom_boxplot( alpha = 0) +
  scale_fill_manual(values= primates_palette) +
  scale_color_manual(values= primates_palette) +
  facet_wrap(~Trait, scales = "free") +
  theme_minimal() + theme(legend.position = "none",
                          axis.title.x  = element_blank(),
                          axis.text.x = element_blank()) + 
  geom_point(aes(x=Trait, 
                 y=value,
                 color = ifelse(outlier != "NA",
                                Family,
                                "white"),
                 shape = ifelse(Species %in% overlapping_species_public,
                                "1",
                                "2")),
             position = position_jitter(height = 0.07,
                                        width = 0.07),
             size = 2.5,
             stroke=1) +
  scale_shape_manual(values = c(8, 16))


#######################################
#######################################
##TREE PLOT FOR PUBLIC SPECIES COLOURED

for (i in 1:length(melted_dataframe_traits$Species)){
  if (melted_dataframe_traits$Species[i] %in% species_names$Cancer_species){
    melted_dataframe_traits$tree_species[i] <- as.character(species_names$Tree_species[which(species_names$Cancer_species == as.character(melted_dataframe_traits$Species[i]))])
  }}

pruned_primates.tree<-drop.tip(tree,tree$tip.label[-match(tree_species[!is.na(tree_species)], 
                                                          tree$tip.label)])

Malignancy_rate <- melted_dataframe_traits$value[melted_dataframe_traits$Trait == "NeoplasiaRate"]
names(Malignancy_rate) <- melted_dataframe_traits$tree_species[melted_dataframe_traits$Trait == "NeoplasiaRate"]

par(fg="transparent")
obj<-contMap(pruned_primates.tree,Malignancy_rate[which(!is.na(names(Malignancy_rate)))],plot=FALSE)
obj<-setMap(obj,c("#07ae00", "#53EC2D", "#F1F97C", 
                  "white", "#F08431", "#F03131"))
lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-lastPP$x.lim
ylim<-lastPP$y.lim
plot(obj,fsize=c(1),lwd=c(5),leg.txt="Neoplasy rate",
     ftype="off", xlim=c(0,1))

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
colors<-ifelse(names(Malignancy_rate) %in% overlapping_species_public,
               "#cc492c",
               "black")
names(colors) <- names(Malignancy_rate) 
for(i in 1:length(colors))
  text(lastPP$xx[i],lastPP$yy[i],obj$tree$tip.label[i],
       pos=4,cex=1,col=colors[obj$tree$tip.label[i]],font=2)



            