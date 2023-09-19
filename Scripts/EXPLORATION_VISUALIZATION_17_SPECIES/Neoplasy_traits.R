#####################
#NEOPLASY RATE#####
#####################

Neoplasia <- cancer_filtered_traits$NeoplasiaRate
names(Neoplasia) <- tree_species
Neoplasia <- Neoplasia[!is.na(names(Neoplasia))]

#Malignant
obj<-contMap(pruned_primates.tree,Neoplasia,plot=FALSE)
obj<-setMap(obj,c("#077223", "#53EC2D", "white", "#F1F97C", 
                  "#F08431", "#F03131"))
plot(obj,fsize=c(1),lwd=c(5),leg.txt="Neoplasia Rate")

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
  mutate(z_score = (value-mean)/sd) %>%
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
  geom_text( fontface = "bold",
             position = position_jitter(height = 0.08,
                                        width = 0.08),
             size=4,
             aes(group = Family, 
                 label = ifelse((label == "high_extreme") |
                                  (label == "low_extreme") ,
                                as.character(Species),'')))  +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white")) +
  ylab("Trait value") +
  guides(fill = guide_legend(override.aes = list(size = 0.5)))

###############
#boxplot trait#
ggplot(data=melted_dataframe_traits, 
       aes(x=Trait, 
           y=value,
           fontface="bold")) +
  geom_boxplot(aes(fill="white"), alpha = 0) +
  facet_wrap(~Trait, scales = "free") +
  #geom_text(aes(label = outlier,
  #            color = Family,
  #           show_legend = FALSE), na.rm = TRUE,
  #      position = position_jitter(height = 0.08,
  #                                width = 0.08)) +
  scale_color_manual(values= primates_palette) +
  scale_fill_manual(values= "white") +
  geom_point(aes(x=Trait, 
                 y=value,
                 color = ifelse(outlier != "NA",
                                Family,
                                "white")),
             position = position_jitter(height = 0.08,
                                        width = 0.08)) +
  theme_minimal()

#boxplot families
ggplot(data=melted_dataframe_traits, 
       aes(x=Trait, 
           y=value,
           color=Family,
           fill=Family,
           fontface="bold")) + geom_boxplot( alpha = 0.6) +
  scale_fill_manual(values= primates_palette) +
  scale_color_manual(values= primates_palette) +
  facet_wrap(~Trait, scales = "free")


 
  

