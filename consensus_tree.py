import os
import subprocess
import re


def run_iqtree(input_dir, filename, output_dir, model='MFP', num_bootstraps=1):
    """
    Runs IQ-TREE with ultrafast bootstrap (ufboot).
    """
    try:
        command = [
            'iqtree',
            '-s', f'{input_dir}/{filename}',
            '-m', model,
            '-b', str(num_bootstraps),  # Ultrafast bootstrap
            '-pre', f'{output_dir}/{filename.strip(".fasta")}',
            '-quiet'
        ]
        subprocess.run(command)

    except subprocess.CalledProcessError as e:
        print(f'Error running IQ-TREE: {e.stderr}')
    
    clean_iqtree_output(output_dir)


def clean_iqtree_output(directory):
    """
    Deletes all files in directory except the ones with .treefile extension.
    """
    command = f'find {directory} -type f ! -name *.treefile -delete'
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f'Error while cleaning up IQ-TREE files: {e}')


def copy_files_to_one(directory, output_file):
    """
    Copies the contents of all files in a directory into a single file.
    """
    try:
        with open(output_file, 'w') as outfile:
            for filename in os.listdir(directory):
                file_path = os.path.join(directory, filename)
                
                if os.path.isfile(file_path):
                    with open(file_path, 'r') as infile:
                        outfile.write(infile.read())
        
        print(f'All files in {directory} have been copied to {output_file}.')
    except Exception as e:
        print(f'Error while copying to file {output_file}: {e}')


def modify_treefile_labels(input_file, output_file):
    """
    Modifies treefile labels by shortening the names to the first four letters.
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                modified_line = re.sub(r'(\b\w+?)(:)', lambda m: m.group(1)[:4] + m.group(2), line)
                outfile.write(modified_line)
        
        print(f'Modified treefile has been saved to {output_file}.')
    except Exception as e:
        print(f'An error occurred: {e}')


def generate_consensus(all_trees_file, consensus_tree_prefix, majority=False):
    """
    Generates a consensus tree from trees in all_trees_file.
    If majority: generates a majority consensus, else greedy consensus.
    """
    try:
        if majority:
            subprocess.run(['iqtree', '-con', all_trees_file, '-pre', consensus_tree_prefix, '-minsup', '0.51'])
        else:
            subprocess.run(['iqtree', '-con', all_trees_file, '-pre', consensus_tree_prefix])
    except Exception as e:
        print(f'Error while generating consensus: {e}')


def generate_super_tree(input_file, output_file):
    """
    Generates a supertree by combining input trees and creating a consensus.
    """
    try:
        subprocess.run(['iqtree', '-stree', input_file, '-pre', output_file])
    except Exception as e:
        print(f'Error while generating supertree: {e}')


if __name__ == '__main__':

    input_dir = 'genolevures_uniquefamilies'
    alignment_dir = 'alignments'
    trees_dir = 'trees'
    all_trees_file = 'all_trees.treefile'
    all_trees_preprocessed_file = 'all_trees_preprocessed.treefile'
    greedy_consensus_prefix = 'greedy_consensus'
    majority_consensus_prefix = 'majority_consensus'
    supertree_prefix = 'supertree'

    for file in os.listdir(input_dir):
        with open(f'{alignment_dir}/{file}', 'w') as alignment_file:
            subprocess.run(['mafft', f'{input_dir}/{file}'], stdout=alignment_file)

    for file in os.listdir(alignment_dir):
        run_iqtree(alignment_dir, file, trees_dir)

    copy_files_to_one(trees_dir, all_trees_file)
    modify_treefile_labels(all_trees_file, all_trees_preprocessed_file)

    generate_consensus(all_trees_preprocessed_file, greedy_consensus_prefix)
    generate_consensus(all_trees_preprocessed_file, majority_consensus_prefix, majority=True)

    generate_super_tree(all_trees_preprocessed_file, supertree_prefix)
