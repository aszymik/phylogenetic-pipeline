import os

input_dir = 'sequences'
output_file = 'blast/proteomes.fasta'

with open(output_file, 'w') as outfile:
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta'):
            filepath = os.path.join(input_dir, filename)
            with open(filepath, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        # Add species name to header
                        outfile.write(f">{filename.strip('.fasta').replace(' ', '_')}_{line[1:]}")
                    else:
                        outfile.write(line)
