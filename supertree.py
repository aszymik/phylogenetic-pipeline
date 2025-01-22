import os
import argparse
import subprocess


def concatenate_trees(input_dir, output_file, extension):
    """
    Concatenates all tree files with the specified extension into one output file.

    Args:
        input_dir (str): Path to the directory containing tree files.
        output_file (str): Path to the output concatenated file.
        extension (str): File extension for tree files (e.g., '.tre', '.nwk').

    Returns:
        None
    """
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(input_dir):
            if filename.endswith(extension):
                filepath = os.path.join(input_dir, filename)
                with open(filepath, 'r') as infile:
                    outfile.write(infile.read())
                    outfile.write('\n')

    print(f"Concatenated all '{extension}' files from {input_dir} into {output_file}.")

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
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing tree files to concatenate.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save the concatenated output file.')
    parser.add_argument('--fasturec_path', type=str, default='fasturec', help='Path to the Fasturec executable. Default: fasturec.')
    parser.add_argument('--extension', type=str, default='.nwk', help='File extension of tree files to concatenate. Default: .nwk')
    args = parser.parse_args()

    # Concatenate tree files
    concatenate_trees(args.input_dir, args.output_file, args.extension)

    # Run Fasturec
    run_fasturec(args.fasturec_path, args.output_file)

if __name__ == '__main__':
    main()
