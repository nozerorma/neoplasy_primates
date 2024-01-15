#!/bin/bash

# Source venv
source /home/miguel/TFM/Master_projects/NEOPLASY_PRIMATES/Scripts/6.Internal_validation/venv/bin/activate

caas_path="5.OncoVar/Results_coincident"
trait_file="Neoplasia_species360/clean_primate_traits.csv"
species_list_path="8.Internal_validation/Other_spp_filtered"
results_path="8.Internal_validation/Results"

mkdir -p ${results_path}

for caas_file in $(find ${caas_path} -name "*.tsv*"); do
    echo $caas_file
    caas_name=$(basename "${caas_file}")
    python Scripts/6.Internal_validation/internal_validation.py ${caas_file} ${trait_file} ${species_list_path}/${caas_name}.filtered_out.tab
done
