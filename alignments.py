import os
import subprocess

def run_mafft_on_families(input_dir, output_dir):
    """
    Aligns sequences in each family using MAFFT.
    """
    for family_file in os.listdir(input_dir):
        output_file = os.path.join(output_dir, f'{family_file[:-6]}_msa.fasta')
        with open(output_file, 'w') as alignment_file:
            subprocess.run(['mafft', f'{input_dir}/{family_file}'], stdout=alignment_file)

if __name__ == '__main__':
    input_dir = 'clustering/clusters'
    output_dir = 'msa'
    os.makedirs(output_dir, exist_ok=True)
    run_mafft_on_families(input_dir, output_dir)
