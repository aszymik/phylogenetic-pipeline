#!/usr/bin/env python3

import os
import subprocess
import argparse


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



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate consenus using greedy consensus and majority consensus from all trees in the given directory.')
    parser.add_argument('--trees_dir', required=True, help='Directory containing input tree files.')
    parser.add_argument('--output_dir', required=True, help='Directory to save generated consensus trees.')
    args = parser.parse_args()

    trees_dir, consensus_dir = args.trees_dir, args.output_dir

    os.makedirs(consensus_dir, exist_ok=True) # create output directory if not exists

    all_trees_file = f'{consensus_dir}/all_trees.treefile'
    greedy_consensus_prefix = f'{consensus_dir}/greedy_consensus'
    majority_consensus_prefix = f'{consensus_dir}/majority_consensus'

    # Copy all trees into one file (required by IQ-TREE)
    copy_files_to_one(trees_dir, all_trees_file)

    # Generate consensus using two methods
    generate_consensus(all_trees_file, greedy_consensus_prefix)
    generate_consensus(all_trees_file, majority_consensus_prefix, majority=True)
    