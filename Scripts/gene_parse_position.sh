#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <file> <species> <position>"
  exit 1
fi

# Assigning command line arguments to variables
file="$1"
species="$2"
position="$3"

# Using awk to parse the file based on species and position
result=$(awk -v species="$species" -v pos="$position" '
$1 == species {print substr($2, pos, 1)}' "$file")

# Checking if the result is empty (species not found or position out of bounds)
if [ -z "$result" ]; then
  echo "Error: Species not found or position out of bounds."
else
  echo "Amino Acid at position $position for $species is: $result"
fi
