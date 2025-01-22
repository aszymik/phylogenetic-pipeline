# Input/output files and directories
bioproject_id = 'PRJNA736147'
sequences_dir = 'sequences'
proteomes_fasta = 'clustering/proteomes.fasta'
mmseqs_results_dir = 'clustering/mmseqs_results'
clusters_dir = 'clustering/clusters'
alignment_dir = 'msa'
trees_dir = 'trees'
consensus_dir = 'consensus'
bootstrap_trees_dir = 'trees_bootstrap'
bootstrap_consensus_dir = 'consensus_bootstrap'

# Parameters used
bootstrap_replicates=100
bootstrap_threshold=70

# Default extensions from scripts
mmseqs_cluster_tsv = f'{mmseqs_results_dir}/mmseqs_cluster.tsv'
all_trees_file = f'{consensus_dir}/all_trees.treefile'
bootstrap_trees_file = f'{bootstrap_consensus_dir}/all_trees.treefile'

# Fetch protein sequences from NCBI database
rule fetch_sequences:
    output:
        directory=sequences_dir
    params:
        bioproject_id=bioproject_id
    script:
        'fetch_sequences.py'

# Merge all sequences into one file
rule merge_into_one_file:
    input:
        directory=sequences_dir
    output:
        fasta=proteomes_fasta
    script:
        'merge_into_one_file.py'

# Cluster sequences using MMSeqs2
rule run_clustering:
    input:
        fasta=proteomes_fasta
    output:
        directory=mmseqs_results_dir
    shell:
        """
        bash clustering.sh {input.fasta} {output.directory}
        """

# Correct clustering to contain orthologs cluster
rule correct_clustering:
    input:
        seq_dir=sequences_dir,
        cluster_tsv=mmseqs_cluster_tsv
    output:
        directory=clusters_dir
    script:
        'correct_clustering.py'

# Run MSA
rule run_alignment:
    input:
        directory=clusters_dir
    output:
        directory=alignment_dir
    script:
        'alignments.py'

# Infer NJ trees
rule infer_trees:
    input:
        directory=alignment_dir
    output:
        directory=trees_dir
    shell:
        """
        Rscript infer_trees.R \
            --alignment_dir {input.directory} \
            --output_dir {output.directory}
        """

# Generate consensus
rule consensus_tree:
    input:
        directory=trees_dir
    output:
        directory=consensus_dir
    script:
        'consensus_tree.py'

# Generate supertree
rule supertree:
    input:
        treefile=all_trees_file
    script:
        'supertree.py'

# Infer NJ trees with bootstrap
rule infer_bootstrap_trees:
    input:
        directory=alignment_dir
    output:
        directory=bootstrap_trees_dir
    params:
        bootstrap_replicates=bootstrap_replicates
        bootstrap_threshold=bootstrap_threshold
    shell:
        """
        Rscript infer_trees.R \
            --alignment_dir {input.directory} \
            --output_dir {output.directory} \
            --use_bootstrap \
            --bootstrap_replicates {params.bootstrap_replicates} \
            --bootstrap_threshold {params.bootstrap_threshold}
        """

# Generate consensus for bootstrapped trees
rule consensus_tree:
    input:
        directory=bootstrap_trees_dir
    output:
        directory=bootstrap_consensus_dir
    script:
        'consensus_tree.py'

# Generate supertree for bootstrapped trees
rule supertree:
    input:
        treefile=bootstrap_trees_file
    script:
        'supertree.py'

