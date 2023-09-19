library(taxize)
library(dplyr)

# Download from ncbi taxonomic rank of defined species list
taxon_list <- tax_name(species_vector, db = "ncbi", accepted=TRUE, get = c("order", "subclass", "family"))

# Check for null values and correct them
any(is.na(taxon_list)) # check NA vals
taxon_list[!complete.cases(taxon_list),] # show NAs

# Correct NAs, fill from rows sharing genus
taxon_list$genus <- sub("^(\\w+)_.*", "\\1", taxon_list$query)# define genus from binomial
unique_genus_values <- unique(taxon_list$genus) # get unique values from genus

taxa_no_na <- function(taxon_list) { #using dplyr
  taxon_list %>%
  group_by(genus) %>%
  mutate(
    order = ifelse(is.na(order), first(order[!is.na(order)]), order),
    suborder = ifelse(is.na(suborder), first(suborder[!is.na(suborder)]), suborder),
    family = ifelse(is.na(family), first(family[!is.na(family)]), family)
  ) %>%
  ungroup()
}

# Create filled taxa df
filled_taxa <- taxa_no_na(taxon_list)

# Create group abbreviation (from suborder)
grouped_taxa <- filled_taxa %>%
  mutate(GROUP = ifelse(suborder == 'Haplorrhini', paste0('HAP_', family), paste0('STR_', family)))

# Save in RDS object
saveRDS(grouped_taxa, "./Data/233-GENOMES/grouped_taxa.rds")

