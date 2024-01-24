# Neoplasy Primates

## Abstract

This repository contains supplementary data and code for the study investigating genetic factors associated with cancer prevalence in primates. The research utilizes phylogenetic comparative methods, specifically CAAS and RER analysis of orthologous proteins. The study uncovers mutations in key cancer driver genes (ABL2, KMT2C, TGFBR2) involved in various carcinogenic pathways. The presence of these mutations across diverse primate species suggests shared evolutionary mechanisms in cancer development. The findings emphasize the importance of analyzing a broader range of genes, considering factors such as body size, longevity, sex, or tumoral origin.

## Repository Content

### Supplementary Data

Supplementary Data may be found under the [General TFM GitHub Repository](https://github.com/nozerorma/neoplasy_primates).

### Supplementary Figures (DATA S1)

1. MSA number per species (figure)
2. Foreground-background CAAS distribution
3. Pvalue CAAS distribution
4. ORA plots
5. Global association assay figures
6. Internal consistency assay figures

### Supplementary Files (DATA S2)

1. R packages citation file
2. Fossil-corrected phylogenetic tree
3. MSA number per species (data)
4. CAAS results
5. RER results
6. CAAS-RER cooccurrence
7. Conda environment configuration
8. Docker container buildfile*

### Scripts
**1.data_preparation.Rmd**

Description: R Markdown file for data preparation.

**2.neoplasy_exploration_def.Rmd**

Description: R Markdown file for neoplasy exploration and analysis.

**3.CAAS_analysis_streamlined.Rmd**

Description: R Markdown file for streamlined CAAS analysis.

**4.Gene_analysis.Rmd**

Description: R Markdown file for gene analysis.

**5.RER_analysis.Rmd**

Description: R Markdown file for RER analysis preview.

**6.Internal_validation**

	Description: R Markdown file for internal consistency analysis (also for association study). Includes: 
 - internal_validation.Rmd
 - internal_validation.py
 - Internal_validation_run.sh

**7.Cross-analysis.Rmd**

Description: R Markdown file for cross-analysis CAAS-RER.

**8.OncoVar_ora_analysis.Rmd**

Description: R Markdown file for OncoVar ORA analysis.

**9.OncoDB_buildup.Rmd**

Description: R Markdown file for OncoDB buildup (extra).

**10.species_apport_alignment**

Description: Python script for analysis of the number of alignments per species. Includes: 
 - species_apport_alignment.r
 - species_apport_alignment.py
 - species_apport_alignment.Rmd

**11.r_package_citation.Rmd**
Description: R Markdown file for citing R packages.

**12.gene_parse_position.sh**
Description: Shell script for parsing gene positions (extra).

## How to Use

To reproduce the analysis or explore the supplementary data, follow the steps below:

1. Clone the repository:

```bash
git clone https://github.com/nozerorma/neoplasy_primates.git
cd neoplasy_primates
```

2. Explore the data, supplementary files and scripts in the respective directories.

3. Use the provided Conda environment configuration and Docker container buildfile for reproducibility and dependency check. 

\* Note that Rocker (Docker) was not used for R analysis. R analysis was performed from the included Conda environment (using a Micromamba installation). Docker installation is included as a remote alternative to the Conda environment, dependencies would need to be installed accordingly.

