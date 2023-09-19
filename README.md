# CANCER PROYECT

## This is a very funky README including the aspects discussed in the first meeting. It wil be updated as soon as things get going.

### ABOUT ACTUAL METHODOLOGY

##### Monovariable assay probably only including BM (Peto's Paradox)

##### Gene Discovery
- **RERConverge**: once per phenotype FOR THE WHOLE PHYLOGENY (nextflow?)
- **CAASTools**: tricky because of the scenarios. Some possibilities that must be taken into account include the double-trouble analysis, first splitting, then co-analyzing (see those nice **Venn Diagrams**)
##### Testing
- **RERC**: no biggies
- **CAASTools**: Bootstrapping **IS** demanding. So do the discovery stage, see if there's something interesting in the in-silico analysis, get to the bootstraping stage.
##### Metadata / In-Silico Funtional Analysis
- GO, KEGG, Pathway enrichment analysis **(TBC)**
- Protein network (**String**)
-- Protein-Protein Interaction Networks
-- Functional Enrichment Analysis
-- Centrality / protein interaction stuff (?)

### ABOUT THE SCENARIOS
Flexibilize this in some way looking at diff combinations (~maybe that sorta eigenvector structures~ through some automatation such as permutations?)
Perfect conservation of AA in top -> two fams (check this out, think it has to do with the restriction criteria, top must be top both family and database-wise)

#### CAAS AND RER
Different scenarios apply to CAAS analysis and will be ~compared~ *added* to the results of RERs.
The basis is doing pretty much what has been done up until now with a more flexible scenario setting leading to a wider gamut of CAAS discoveries (not for sure)
