### A very small script to prune alignment apport using only our dataset species

library(phytools)
library(dplyr)
library(palettetown)

ali_apport <- read.csv("../Data/alignment_apport.tsv", sep = "\t")

cancer_traits <- read.csv("../Data/Neoplasia_species360/clean_primate_traits.csv")

tree <- read.tree("../Data/Phylo-233-GENOMES/science.abn7829_data_s4.nex.tree")
tree_species <- tree$tip.label


common_names <- intersect(tree$tip.label, cancer_traits$species)

### Mantain structure

### Mantain structure
ali_apport <- ali_apport %>%
  mutate(ordering = match(Pattern, common_names)) %>%
  arrange(ordering)  %>%
  filter(!is.na(ordering)) %>%
  select(!ordering) %>%
  rename(Species = Pattern) %>%
  select(Family, everything()) %>%
  arrange(Family)

write.csv(ali_apport, "../Data/ali_apport.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


ordered_species <- ali_apport$Species
ali_apport$Species <- factor(ali_apport$Species, levels = ordered_species)

# Assuming your data is stored in a data frame called 'ali_apport'
ggplot(ali_apport, aes(x = Family, y = Count, fill = Family)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Counts per Species Grouped by Family",
       x = "Species",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


scatter_plot <- ggplot(data = ali_apport, 
                       aes(x = Species, y = Count, color = Family)) +
  geom_point(size = 3) +
  geom_line(aes(group = Family), size = 1, alpha = 0.5) +
  labs(x = "Primate Species", y = "Number of alignments") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 0, angle = 90, hjust = 1),
    axis.title.x = element_text(size = 16, hjust = 0.5, margin = margin(t = 20, b = 20)),
    axis.title.y = element_text(size = 16, angle = 90, hjust = 0.5, margin = margin(l = 20, r = 20)),
    plot.title = element_text(size = 16, face = "bold", color = pokepal("rayquaza", spread = 1), hjust = 0.5, margin = margin(t = 20, b = 20)),
    plot.caption = element_text(size = 10, face = "bold", margin = margin(t = 20, b = 20)),
    legend.position = "right",
    legend.key.size = unit(1, "cm"), legend.text = element_text(size = 15), legend.title = element_text(size = 15))

scatter_plot

ggsave("../Data/ali_apport_scatter.png", scatter_plot, height = 5, width = 10, dpi = 300, dev = "png")

# Convert Species to a factor with a specific order based on Family

library(ggplot2)

# Convert Species to a factor with a specific order based on Family
