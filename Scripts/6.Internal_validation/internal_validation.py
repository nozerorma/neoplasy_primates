from Bio import AlignIO
import pandas as pd
import os
import re
import glob
import sys
import pdb


# Define your data directories
dataDir = '/home/miguel/IBE-UPF/TFM/Master_projects/NEOPLASY_PRIMATES/Data'
outDir = '/home/miguel/IBE-UPF/TFM/Master_projects/NEOPLASY_PRIMATES/Out'

# Expand the user's home directory in the directories path
expanded_dataDir = os.path.expanduser(dataDir)
expanded_outDir = os.path.expanduser(outDir)

# Check if the file path argument is provided
if len(sys.argv) < 4:
    print("Please provide a CAAS file path, a traitfile path and species file path as arguments.")
    sys.exit(1)

# Assign the relative file path from the first command-line argument
caasfile_relative_file_path = sys.argv[1]
traitfile_relative_file_path = sys.argv[2]

# File paths for the dataframes saved from R
# File for families contrast
# Combine the base directory with the relative file path
caasfile_path = os.path.join(expanded_outDir, caasfile_relative_file_path)

traitfile_path = os.path.join(expanded_dataDir, traitfile_relative_file_path)

# Read the dataframe into a pandas DataFrame
df_byFGN = pd.read_csv(caasfile_path, sep="\t")

traits = pd.read_csv(traitfile_path)
# Remove rows where 'neoplasia_prevalence' column has NaN or missing values
traits = traits.dropna(subset=['neoplasia_prevalence'])

# Function to process substitutions with multiple amino acids
def long_substitution(substitution):
    # Extract all letters before the slash using regex
    long_aa = re.findall(r"([A-Z]+)/", substitution)[0]
    return list(long_aa) if len(long_aa) > 1 else [long_aa]

def short_substitution(substitution):
    # Extract all letters after the slash using regex
    matches = re.findall(r"[A-Z]+/([A-Z]+)", substitution)
    if not matches:
        # Handle the case where there is no match
        return None
    short_aa = matches[0]  # Get the first (and likely only) match
    return list(short_aa) if len(short_aa) > 1 else [short_aa]

# Define the species and the substitutions of interest

# Assign the other species file path from the command-line argument
other_spp_filepath = sys.argv[3]
other_spp_file_path = os.path.join(expanded_outDir, other_spp_filepath)

# Initialize an empty list to store the species names
intermediate_quantile_spp = []

# Open the file and read each line
with open(other_spp_file_path, 'r') as file:
    for line in file:
        # Add the species name to the list, stripping any trailing newline characters
        intermediate_quantile_spp.append(line.strip())

# Create a list to hold the data for the new DataFrame
species_classification_data = []


# Check alignment by alignment 
# Only iterating over multiple alignments from genes in which a CAAS was found a CAAS
for gene in df_byFGN['Gene_symbol'].unique():
    # Extract substitutions for the current gene
    gene_data = df_byFGN[df_byFGN['Gene_symbol'] == gene]
    
    # Initialize the substitutions dictionaries
    substitutions_long = {}
    substitutions_short = {}

    for idx, row in gene_data.iterrows():
        # If the position is already in the dictionary, append the amino acid(s)
        position = row['discovery.Position']

        long_aas = long_substitution(row['discovery.Substitution'])
        if position in substitutions_long:
            substitutions_long[position] += long_aas
        else:
            substitutions_long[position] = long_aas

        short_aas = short_substitution(row['discovery.Substitution'])
        if position in substitutions_short:
            substitutions_short[position] += short_aas
        else:
            substitutions_short[position] = short_aas
            
    # Read the alignment file
    pattern = os.path.join(expanded_dataDir, f"Primate_alignments/{gene}.*.phy")

    # Use glob to find files that match the pattern
    matching_files = glob.glob(pattern)

    # Now you can iterate over all matching files
    for file_path in matching_files:
        try:
            alignment = AlignIO.read(file_path, 'phylip-relaxed')
            # Process the alignment as before
        except Exception as e:
            print(f"An error occurred while processing the file {file_path}: {e}")

    # Check for each substitution in the species of interest
    for record in alignment:
        if record.id in intermediate_quantile_spp:
            print(intermediate_quantile_spp)
            for position in gene_data['discovery.Position'].unique():
                # Default classification
                trait_classification = None
                neoplasia_prevalence = None

                # Check long substitutions
                if position in substitutions_long and record.seq[position] in substitutions_long[position]:
                    print(f"CAAS {position} {gene_data.loc[gene_data['discovery.Position'] == position, 'discovery.Substitution'].values[0]} in {gene}, {record.id} has low-NP related AA")
                    trait_classification = 'Low-NP'
                    neoplasia_prevalence = traits.loc[traits['species'] == record.id, 'neoplasia_prevalence'].values[0]

                # Check short substitutions
                if position in substitutions_short and record.seq[position] in substitutions_short[position]:
                    print(f"CAAS {position} {gene_data.loc[gene_data['discovery.Position'] == position, 'discovery.Substitution'].values[0]} in {gene}, {record.id} has high-NP related AA")
                    trait_classification = 'High-NP'
                    neoplasia_prevalence = traits.loc[traits['species'] == record.id, 'neoplasia_prevalence'].values[0]

                # Ensure we only add records where a classification was made
                if trait_classification:
                    # Add the data to our list
                    species_classification_data.append({
                        'Gene': gene,
                        'Position': position,
                        'Substitution': gene_data.loc[gene_data['discovery.Position'] == position, 'discovery.Substitution'].values[0],
                        'label': record.id,
                        'type_neoplasia_prevalence': trait_classification,
                        'neoplasia_prevalence': neoplasia_prevalence
                        })

# Create a DataFrame from the collected data
df_classification = pd.DataFrame(species_classification_data)

# Get the prefix from the first command-line argument
file_name_parts = os.path.basename(sys.argv[1]).split('_')
if len(file_name_parts) >= 2:
    file_prefix = '_'.join(file_name_parts[:2])
else:
    # Handle case where there's only one part or none
    file_prefix = file_name_parts[0] if file_name_parts else "default"

# Append the prefix to the output file name
output_file_name = f"{file_prefix}_sppclassification_results.tab"

# Define the path for the output CSV file
output_file_path = os.path.join(expanded_outDir, f"8.Internal_validation/Oncovar/{output_file_name}")

# Export the DataFrame to a tab-separated file
df_classification.to_csv(output_file_path, sep='\t', index=False)

print(f"Results have been saved to {output_file_path}")

