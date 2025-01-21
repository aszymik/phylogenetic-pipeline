import os
import subprocess

def run_mafft_on_families(input_dir, output_dir):
    """
    Aligns sequences in each family using MAFFT.
    """
    for folder in os.listdir(input_dir):
        family_dir = os.path.join(input_dir, folder)
        if os.path.isdir(family_dir):
            output_file = os.path.join(output_dir, f'{folder}_aligned.fasta')
            with open(output_file, 'w') as alignment_file:
                subprocess.run(['mafft', '--auto', '--reorder', family_dir], stdout=alignment_file)
            print(f'Alignment for family {folder} completed.')

if __name__ == '__main__':
    input_dir = 'clustering/clusters'
    output_dir = 'multi_alignments'
    os.makedirs(output_dir, exist_ok=True)
    run_mafft_on_families(input_dir, output_dir)
