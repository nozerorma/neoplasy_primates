import os
import re

def count_occurrences_in_files(directory):
    occurrences_count = {}

    # Compile regular expression pattern for "Name1_name2"
    pattern = re.compile(r'\b([A-Za-z0-9]+_[A-Za-z0-9]+)\b')

    # Iterate through all files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith(".phy"):
            file_path = os.path.join(directory, filename)

            # Open and read the content of the file
            with open(file_path, 'r') as file:
                file_content = file.read()

                # Find all occurrences of the pattern using regular expression
                matches = re.findall(pattern, file_content)

                # Update the occurrences count for each match
                for match in matches:
                    occurrences_count[match] = occurrences_count.get(match, 0) + 1

    return occurrences_count

def read_species_family_correspondence(correspondence_file):
    species_to_family = {}

    with open(correspondence_file, 'r') as file:
        for line in file:
            species, family = line.strip().split('\t')
            species_to_family[species] = family

    return species_to_family

def write_occurrences_with_family_to_tsv(occurrences_count, species_to_family, output_file):
    with open(output_file, 'w') as tsv_file:
        # Write header
        tsv_file.write("Pattern\tCount\tFamily\n")

        # Write data
        for pattern, count in occurrences_count.items():
            family = species_to_family.get(pattern, 'Unknown')
            tsv_file.write(f"{pattern}\t{count}\t{family}\n")


# Example usage:
directory_path = '/home/miguel/IBE-UPF/TFM/Master_projects/NEOPLASY_PRIMATES/Data/Primate_alignments'
output_tsv_file = '/home/miguel/IBE-UPF/TFM/Master_projects/NEOPLASY_PRIMATES/Data/alignment_apport.tsv'
correspondence_file = '/home/miguel/IBE-UPF/TFM/Master_projects/NEOPLASY_PRIMATES/Data/Neoplasia_species360/sp2fam_wo_groups.tab'

occurrences_count = count_occurrences_in_files(directory_path)

# Print the total occurrences for each pattern
for pattern, count in occurrences_count.items():
    print(f"Total occurrences of '{pattern}': {count}")

# Read species to family correspondence
species_to_family = read_species_family_correspondence(correspondence_file)

# Write occurrences with family information to TSV file
write_occurrences_with_family_to_tsv(occurrences_count, species_to_family, output_tsv_file)

print(f"Total occurrences with family information written to {output_tsv_file}")
