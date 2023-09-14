######################################
####All species trait distribution####
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

primates_palette <- readRDS("~/Documents/PhD_EvoGenomics/1st_year/GenomePhenome/src/primate_families_palette.rds")



###############################################
##Now let's see what happens with the traits###
###############################################

cancer_species_traits <- read.csv("../data/Phenome/cancer_traits.csv", header=T, sep = "\t")

hist(cancer_species_traits$PropMalignant)
#total variance analysis

for (i in 1:length(cancer_species_traits$Species)){
  if(as.character(cancer_species_traits$Species[i]) %in% all_traits$SpeciesBROAD){
    cancer_species_traits$Group[i] <- as.character(all_traits$groups[as.character(all_traits$SpeciesBROAD) == as.character(cancer_species_traits$Species[i])])
  }
  else{
    cancer_species_traits$Group[i] <- NA
  }}

#MELT dataset

library(e1071)
skewness(cancer_species_traits$PropMalignant)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}





#TREE WITH INFORMATION

#tree
tree <- read.newick("../data/Tree/primates.newick")
#species names
species_names <- read.csv("../data/Phenome/names.txt", header=T, sep = "\t")


tree_species <- c()

###################################
#ANALYSIS WITH ONLY TREE SPECIES#
###################################

cancer_filtered_traits <- cancer_species_traits %>% 
  group_by(Species) %>%
  summarise_all(funs(mean(., na.rm=TRUE)))
View(cancer_filtered_traits)

tree_species <- c()
for (i in 1:length(cancer_filtered_traits$Species)){
  if (cancer_filtered_traits$Species[i] %in% species_names$Cancer_species){
    tree_species[i] <- as.character(species_names$Tree_species[which(species_names$Cancer_species == as.character(cancer_filtered_traits$Species[i]))])
  }}

pruned_primates.tree<-drop.tip(tree,tree$tip.label[-match(tree_species[!is.na(tree_species)], 
                                                          tree$tip.label)])

#####################
#MALIGNANCY RATE#####
#####################

Malignancy <- cancer_filtered_traits$MalignancyRate
names(Malignancy) <- tree_species
Malignancy <- Malignancy[!is.na(names(Malignancy))]

#Malignant
obj<-contMap(pruned_primates.tree,Malignancy,plot=FALSE)
obj<-setMap(obj,c("#077223", "#53EC2D", "white", "#F1F97C", 
                  "#F08431", "#F03131"))
plot(obj,fsize=c(1),lwd=c(5),leg.txt="Malignancy Rate")


#boxplot
for (i in 1:length(cancer_filtered_traits$Species)){
  if(as.character(cancer_filtered_traits$Species[i]) %in% all_traits$SpeciesBROAD){
    cancer_filtered_traits$Group[i] <- as.character(all_traits$groups[as.character(all_traits$SpeciesBROAD) == as.character(cancer_filtered_traits$Species[i])])
  }
  else{
    cancer_filtered_traits$Group[i] <- NA
  }}


medians <- t(cancer_filtered_traits[,-which(colnames(cancer_filtered_traits) %in% c("Species", "Group"))] %>% 
  summarize_all(., median))

means <- t(cancer_filtered_traits[,-which(colnames(cancer_filtered_traits) %in% c("Species", "Group"))] %>% 
             summarize_all(., mean))

sds <- t(cancer_filtered_traits[,-which(colnames(cancer_filtered_traits) %in% c("Species", "Group"))] %>% 
          summarize_all(., funs(sd)))

maxs <-  t(cancer_filtered_traits[,-which(colnames(cancer_filtered_traits) %in% c("Species", "Group"))] %>% 
                   summarize_all(., max))

melted_dataframe_traits <- melt(cancer_filtered_traits, id.vars = c("Group", "Species"))
colnames(melted_dataframe_traits) <- c("Family", "Species", "Trait", "value")

for (i in 1:length(melted_dataframe_traits$Species)){
  if(as.character(melted_dataframe_traits$Trait[i]) %in% rownames(medians)){
    melted_dataframe_traits$median[i] <- as.numeric(as.character(medians[rownames(medians) == melted_dataframe_traits$Trait[i]]))
    melted_dataframe_traits$mean[i] <- as.numeric(as.character(means[rownames(means) == melted_dataframe_traits$Trait[i]]))
    melted_dataframe_traits$sd[i] <- as.numeric(as.character(sds[rownames(sds) == melted_dataframe_traits$Trait[i]]))
    melted_dataframe_traits$max[i] <- as.numeric(as.character(maxs[rownames(maxs) == melted_dataframe_traits$Trait[i]]))
    
     }}


#all together
melted_dataframe_traits <- melted_dataframe_traits %>%
  dplyr::group_by(Trait) %>% 
  mutate(label = case_when(
    value < quantile(value, 0.25, na.rm=TRUE) ~ "low_extreme", 
      value > quantile(value, 0.75, na.rm=TRUE) ~ "high_extreme",
    TRUE ~ "normal")
  )

View(melted_dataframe_traits)


#x family
melted_dataframe_traits <- melted_dataframe_traits %>%
dplyr::group_by(Family, Trait) %>% 
  mutate(z_score = (value-mean)/sd) %>%
  mutate(label = case_when(
    value < 0.25*max | 
      value > max-(0.25*max) ~ "extreme",
    TRUE ~ "normal")
  )

ggplot(data=melted_dataframe_traits[which(melted_dataframe_traits$Trait == "PropMalignant"),], 
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
  geom_text( fontface = "bold",
             position = position_jitter(height = 0.05,
                                        width = 0.05),
    size=4,
         aes(group = Family, 
              label = ifelse(label == "high_extreme",
                            as.character(Species),'')))  +
theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white")) +
  ylab("Trait value") +
  guides(fill = guide_legend(override.aes = list(size = 0.5)))

#######################
#CORRELATIONS ANALYSIS#
#######################

var_quanti <- PCA(log(cancer_filtered_traits[,2:17]+0.001), 
                  scale.unit = TRUE)
fviz_pca_var(var_quanti, col.var="contrib",
             repel = TRUE)
var_quanti_cor <- as.data.frame(var_quanti$var$cor)
library(Hmisc)
library(corrplot)
res <- cor(log(cancer_filtered_traits[,2:17]+0.001),
           use="pairwise.complete.obs")
corrplot(res,
         tl.col = "black")


































#Coefficient of variation analysis
GrouppedCancer_traits_df <- melt(cancer_filtered_traits, id.vars = c("Species",
                                                                  "Group"))
ContinousTraitsMeans_clades <- GrouppedCancer_traits_df %>%
  group_by(Group, variable) %>%
  summarize(mean(value, na.rm = TRUE))

ContinousTraitsSds_clades <- GrouppedCancer_traits_df %>%
  group_by(Group, variable) %>%
  summarize(sd(value, na.rm = TRUE))

sd(cancer_filtered_traits$PropMalignant[which(cancer_filtered_traits$Group == "STR_Indriidae")])

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
