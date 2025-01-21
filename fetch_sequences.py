#!/usr/bin/env python3

from Bio import Entrez
import os
import argparse

Entrez.email = 'email@example.com'


def get_proteomes(bioproject_id, output_dir):
    """
    Fetch protein sequences from a given BioProject ID and save them to the output directory.

    Args:
        bioproject_id (str): BioProject ID to search for protein sequences.
        output_dir (str): Directory to save the fetched protein FASTA files.
    """
    # Search for assemblies in the given BioProject
    search_handle = Entrez.esearch(db='assembly', term=bioproject_id, retmax=1000)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    assembly_ids = search_results['IdList']
    print(f'Found {len(assembly_ids)} assemblies in BioProject {bioproject_id}.')

    for assembly_id in assembly_ids:
        try:
            # Fetch assembly record to find protein links
            summary_handle = Entrez.esummary(db='assembly', id=assembly_id)
            summary = Entrez.read(summary_handle)
            summary_handle.close()

            species = summary['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName'].replace(' ', '_')

            # FTP link for protein sequences
            ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
            if not ftp_path:
                print(f'No GenBank data available for assembly {assembly_id} ({species}).')
                continue

            protein_fasta_url = f'{ftp_path}/{ftp_path.split("/")[-1]}_protein.faa.gz'
            print(f'Downloading proteins from {protein_fasta_url}...')

            # Download and unzip fasta
            os.system(f'wget -q {protein_fasta_url} -O {output_dir}/{species}.faa.gz')
            os.system(f'gunzip -f {output_dir}/{species}.faa.gz')

            print(f'Downloaded sequences for assembly {assembly_id} ({species}).')

        except Exception as e:
            print(f'An error occurred for assembly {assembly_id}: {e}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch protein sequences from a BioProject ID.")
    parser.add_argument("--bioproject_id", required=True, help="BioProject ID to fetch sequences from.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the protein FASTA files.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)  # ensure output directory exists
    get_proteomes(args.bioproject_id, args.output_dir)  # fetch proteomes
