import pandas as pd
import subprocess

blast_results_path = 'blast/blast_results.tsv'
blast_graph_path = 'blast/blast_graph.mci'
mcl_results_path = 'clustering/clusters.txt'

# Read BLAST data
blast_results = pd.read_csv(blast_results_path, sep='\t', header=None)
blast_results.columns = ['QueryID', 'SubjectID', 'Identity', 'Alignment_length', 'Mismatches', 'Gaps', 
                         'Algn_start_q', 'Algn_end_q', 'Algn_start_s', 'Algn_end_s', 'E_value', 'Bit_score']

# Normalize bit score results to [0, 1]
blast_results['Weight'] = ((blast_results['Bit_score'] - blast_results['Bit_score'].min()) 
                           / (blast_results['Bit_score'].max() - blast_results['Bit_score'].min()))

# Save to MCL format file
blast_results[['QueryID', 'SubjectID', 'Weight']].to_csv(
    'blast/blast_graph.mci', sep=' ', header=False, index=False
)

# Run MCL
# I: inflation parameter
subprocess.run(['mcl', blast_graph_path, '--abc', '-I', '4.0', '-o', mcl_results_path])
