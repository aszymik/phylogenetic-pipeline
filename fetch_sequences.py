from Bio import Entrez, SeqIO
import ssl
import os

Entrez.email = 'email@example.com'
ssl._create_default_https_context = ssl._create_unverified_context

def get_protein_sequences(species_name):
    try:
        # Step 1: Search for the complete genome
        search_term = f'SARS {species_name} complete genome'
        handle = Entrez.esearch(db='nucleotide', term=search_term, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()
        
        if not search_results['IdList']:
            print(f'No genome found for {species_name}')
            return None
        
        genome_id = search_results['IdList'][0]
        print(f'Found genome ID: {genome_id}')
        
        # Step 2: Fetch the GenBank record
        handle = Entrez.efetch(db='nucleotide', id=genome_id, rettype='gb', retmode='text')
        record = SeqIO.read(handle, 'genbank')
        handle.close()
        
        # Step 3: Extract protein ids and sequences
        proteins = {}
        for feature in record.features:
            if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                protein_id = feature.qualifiers['protein_id'][0]
                protein_seq = feature.qualifiers['translation'][0]
                proteins[protein_id] = protein_seq
        
        print(f'Extracted {len(proteins)} protein sequences for {species_name}')
        return proteins

    except Exception as e:
        print(f'An error occurred: {e}')
        return None



with open('species_ids.txt') as species_file: 
    species_names = [name.strip() for name in species_file.readlines()]   


for species_name in species_names:
    protein_sequences = get_protein_sequences(species_name)

    # Save to a FASTA file (optional)
    if protein_sequences:
        with open(f'sequences/{species_name}.fasta', 'w') as fasta_file:
            for protein_id, sequence in protein_sequences.items():
                fasta_file.write(f'>{protein_id}\n{sequence}\n')
        print(f'Protein sequences saved to {species_name}_proteins.fasta')




# def fetch_proteome(species_name, output_dir):
#     '''Fetch the proteome for a given species from the NCBI database.'''
#     try:
#         # Search for the proteome
#         search_term = f'SARS {species_name} complete genome'
#         handle = Entrez.esearch(db='nucleotide', term=search_term, retmax=1)
#         record = Entrez.read(handle)
#         handle.close()

#         if not record['IdList']:
#             print(f'No proteome found for {species_name}.')
#             return

#         # Fetch the protein sequences
#         protein_id = record['IdList'][0]
#         handle = Entrez.efetch(db='protein', id=protein_id, rettype='fasta', retmode='text')
#         proteome_data = handle.read()
#         handle.close()

#         # Save to a file
#         output_file = os.path.join(output_dir, f'{species_name.replace(' ', '_')}.fasta')
#         with open(output_file, 'w') as f:
#             f.write(proteome_data)
#         print(f'Proteome for {species_name} saved to {output_file}.')

#     except Exception as e:
#         print(f'An error occurred while fetching the proteome for {species_name}: {e}')


# def main(input_file, output_dir):
#     '''Fetch proteomes for all species listed in the input file.'''
#     # Create output directory if it doesn't exist
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     # Read species names from the input file
#     with open(input_file, 'r') as f:
#         species_names = [line.strip() for line in f if line.strip()]

#     # Fetch proteomes for each species
#     for species_name in species_names:
#         print(f'Fetching proteome for {species_name}...')
#         fetch_proteome(species_name, output_dir)


# if __name__ == '__main__':
#     input_file = 'species_names.txt'  # Replace with your input file name
#     output_dir = 'proteomes'  # Replace with your desired output directory
#     main(input_file, output_dir)
