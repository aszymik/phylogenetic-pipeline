#!/usr/bin/env python3

import os
import subprocess
import argparse

def run_mafft_on_families(input_dir, output_dir):
    """
    Aligns sequences in each cluster/family using MAFFT.

    Args:
        input_dir (str): Directory containing cluster FASTA files.
        output_dir (str): Directory to store the resulting MSAs.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    for family_file in os.listdir(input_dir):
        # Only process .fasta files or any extension that fits your clusters
        if family_file.endswith('.fasta'):
            output_file = os.path.join(output_dir, f'{family_file[:-6]}_msa.fasta')
            with open(output_file, 'w') as alignment_file:
                subprocess.run(['mafft', f'{os.path.join(input_dir, family_file)}'], 
                               stdout=alignment_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run MAFFT on cluster/family files.')
    parser.add_argument('--input_dir', required=True, help='Directory containing cluster FASTA files.')
    parser.add_argument('--output_dir', required=True, help='Directory to store the resulting MSAs.')
    args = parser.parse_args()

    run_mafft_on_families(args.input_dir, args.output_dir)
