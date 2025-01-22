import os
import re
import argparse
import subprocess


def clean_newick(input_file, output_file):
    """
    Removes invalid negative branch lengths and other parsing issues from a Newick file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Remove negative branch lengths (e.g., '-0.009203788108' -> '0.0')
            cleaned_line = re.sub(r'-\d+(\.\d+)?', '0.0', line)
            outfile.write(cleaned_line)
    print(f"Cleaned Newick file saved to: {output_file}")


def run_fasturec(fasturec_path, input_file):
    """
    Runs Fasturec with the specified input file and the -Z option.

    Args:
        fasturec_path (str): Path to the Fasturec executable.
        input_file (str): Path to the concatenated input file for Fasturec.

    Returns:
        None
    """
    try:
        subprocess.run([fasturec_path, '-G', input_file, '-Z'], check=True)
        print(f'Fasturec successfully ran on {input_file}.')
    except subprocess.CalledProcessError as e:
        print(f'Error while running Fasturec: {e}')

def main():
    parser = argparse.ArgumentParser(description='Concatenate tree files and generate a supertree using Fasturec program.')
    parser.add_argument('--trees_file', type=str, required=True, help='File containing tree files to concatenate.')
    parser.add_argument('--fasturec_path', type=str, default='fasturec', help='Path to the Fasturec executable. Default: fasturec.')
    parser.add_argument('--extension', type=str, default='.nwk', help='File extension of tree files to concatenate. Default: .nwk')
    args = parser.parse_args()

    all_trees_file = args.trees_file
    all_trees_file_cleaned = f'{all_trees_file}_cleaned.nwk'

    # Remove negative edges from file
    clean_newick(all_trees_file, all_trees_file_cleaned)

    # Run Fasturec
    run_fasturec(args.fasturec_path, all_trees_file_cleaned)

    
if __name__ == '__main__':
    main()
