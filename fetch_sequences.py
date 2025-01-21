from Bio import Entrez
import os

Entrez.email = 'email@example.com'


def get_proteoms(bioproject_id, output_dir):
    """Fetch protein sequences from given BioProject"""

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

            protein_fasta_url = f'{ftp_path}/{ftp_path.split('/')[-1]}_protein.faa.gz'
            print(f'Downloading proteins from {protein_fasta_url}...')

            # Download and unzip fasta
            os.system(f'wget -q {protein_fasta_url} -O {output_dir}/{species}.faa.gz')
            os.system(f'gunzip -f {output_dir}/{species}.faa.gz')
        
            print(f'Downloaded sequences for assembly {assembly_id} ({species}).')
        
        except Exception as e:
            print(f'An error occurred for assembly {assembly_id} ({species}): {e}')


output_dir = 'drosophila_sequences'
bioproject_id = 'PRJNA736147'

os.makedirs(output_dir, exist_ok=True)
get_proteoms(bioproject_id, output_dir)
