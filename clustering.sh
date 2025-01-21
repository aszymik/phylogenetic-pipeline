#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta> <output_dir>"
    exit 1
fi

# Input arguments
input_fasta=$1
output_dir=$2

# Create output and temporary directories if they don't exist
mkdir -p "${output_dir}/mmseqs"
mkdir -p "${output_dir}/tmp"

# Run the mmseqs easy-linclust command
mmseqs easy-linclust "$input_fasta" "${output_dir}/mmseqs" "${output_dir}/tmp" --min-seq-id 0.35

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "MMseqs clustering completed successfully."
else
    echo "MMseqs clustering failed."
    exit 1
fi
