import os
import re

input_dir = 'sequences'
output_file = 'blast/proteomes.fasta'


def transform_description(description):
    # Extract the species name within square brackets
    species_name_match = re.search(r'\[(.*?)\]', description)
    if not species_name_match:
        raise ValueError('Species name not found in square brackets.')
    species_name = species_name_match.group(1).replace(' ', '_')
    
    # Extract the protein ID at the beginning
    protein_id_match = re.match(r'>(\S+)', description)
    if not protein_id_match:
        raise ValueError('Protein ID not found at the beginning of the string.')
    protein_id = protein_id_match.group(1)
    
    return f'>{species_name}:{protein_id}\n'

with open(output_file, 'w') as outfile:
    for filename in os.listdir(input_dir):
        # if filename.endswith('.fasta'):
        if filename.endswith('.faa'):
            filepath = os.path.join(input_dir, filename)
            with open(filepath, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        # Add species name to header
                        outfile.write(transform_description(line))
                    else:
                        outfile.write(line)
