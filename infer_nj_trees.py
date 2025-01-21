import subprocess
import os

def run_mafft_alignment(input_dir, output_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta'):
            filepath = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, f'{filename[:-6]}_msa.fasta')

        with open(output_path, 'w') as out_file:
            subprocess.run(['mafft', filepath], stdout=out_file)
            
        subprocess.run(['Rscript', 'infer_trees.R'])