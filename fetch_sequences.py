from Bio import Entrez, SeqIO
import ssl
import os

Entrez.email = 'email@example.com'
ssl._create_default_https_context = ssl._create_unverified_context

no_genome = []
no_proteome = []

def get_protein_sequences(species_name):
    try:
        # Step 1: Search for the complete genome
        search_term = f'{species_name} complete genome'
        handle = Entrez.esearch(db='nucleotide', term=search_term, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()
        
        if not search_results['IdList']:
            print(f'No genome found for {species_name}')
            no_genome.append(species_name)
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
        print(f'Protein sequences saved to {species_name}.fasta')
    else:
        no_proteome.append(species_name)

print(f'NO GENOME:\n{no_genome}')
# 'YN2011'

print(f'NO PROTEOME:\n{no_proteome}')
# ['SARS SZ3', 'PC4-227', 'YN2011']