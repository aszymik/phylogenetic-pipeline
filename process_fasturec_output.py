import argparse


def process_fasturec_output(input_file, output_file):
    """
    Processes a FASTUREC output file by removing the leading number and whitespace,
    and appending a semicolon at the end of the line.

    Args:
        input_file (str): Path to the input file containing FASTUREC output.
        output_file (str): Path to the output file to save the processed data.
    """
    try:


        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                # Remove leading number and whitespace
                processed_line = line.lstrip().split(' ', 1)[-1].strip()
                # Add semicolon at the end if not already present
                if not processed_line.endswith(';'):
                    processed_line += ';'
                # Write the processed line to the output file
                outfile.write(processed_line + '\n')

        print(f"Processed FASTUREC output saved to: {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Processes a FASTUREC output file by removing the leading number and whitespace, and appending a semicolon at the end of the line.')
    parser.add_argument('--fasturec_tree', required=True, help='Path to a file containing Fasturec output tree you want to process.')
    parser.add_argument('--output_tree', required=True, help='Desired output file path.')
    args = parser.parse_args()

    process_fasturec_output(args.fasturec_tree, args.output_tree)
