#!/usr/bin/env python3

import os
import re
import argparse


def transform_description(description):
    """
    Transforms a FASTA header from NCBI-fetched files.

    Args:
        description (str): The header line of a FASTA entry.

    Returns:
        str: Modified header line in the format '>species_name:protein_id'.
    """
    # Extract the species name from square brackets
    species_name_match = re.search(r'\[(.*?)\]', description)
    if not species_name_match:
        raise ValueError('Species name not found in square brackets.')
    species_name = species_name_match.group(1).replace(' ', '_')

    # Extract the protein ID at the beginning of header
    protein_id_match = re.match(r'>(\S+)', description)
    if not protein_id_match:
        raise ValueError('Protein ID not found at the beginning of the string.')
    protein_id = protein_id_match.group(1)

    return f'>{species_name}:{protein_id}\n'


def merge_fasta_files(input_dir, output_file):
    """
    Merges all FASTA files in the input directory into a single output file.

    Args:
        input_dir (str): Directory containing input FASTA files.
        output_file (str): Path to the merged output FASTA file.
    """
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(input_dir):
            if filename.endswith('.faa'):  # .faa extension from NCBI
                filepath = os.path.join(input_dir, filename)
                with open(filepath, 'r') as infile:
                    for line in infile:
                        if line.startswith('>'):
                            # Transform the header line
                            outfile.write(transform_description(line))
                        else:
                            # Write sequence lines as is
                            outfile.write(line)
    print(f"Merged FASTA files into {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge all FASTA files in a directory into one file.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTA files.")
    parser.add_argument("--output_file", required=True, help="Path to the merged output FASTA file.")
    args = parser.parse_args()

    merge_fasta_files(args.input_dir, args.output_file)
