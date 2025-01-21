#!/bin/bash
# Script to perform clustering using MMseqs

input_fasta=$1
output_dir=$2

# Create directories if they don't exist
mkdir -p "${output_dir}/mmseqs"
mkdir -p "${output_dir}/tmp"

# Run the MMseqs easy-linclust command
mmseqs easy-linclust "$input_fasta" "${output_dir}/mmseqs" "${output_dir}/tmp" --min-seq-id 0.35

