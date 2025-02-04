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

clusters_paralogs_dir = 'clustering/clusters_paralogs'
paralogs_alignment_dir = 'msa_paralogs'
paralogs_trees_dir = 'trees_paralogs'
paralogs_consensus_dir = 'consensus_paralogs'

# Parameters used
bootstrap_replicates=100
bootstrap_threshold=70

# Default extensions from scripts
mmseqs_cluster_tsv = f'{mmseqs_results_dir}/mmseqs_cluster.tsv'
all_trees_file = f'{consensus_dir}/all_trees.treefile'
bootstrap_trees_file = f'{bootstrap_consensus_dir}/all_trees.treefile'
paralogs_trees_file = f'{paralogs_consensus_dir}/all_trees.treefile'


rule all:
    input:
        sequences_dir,
        proteomes_fasta,
        mmseqs_results_dir,
        mmseqs_cluster_tsv,
        clusters_dir,
        alignment_dir,
        trees_dir,
        consensus_dir,


# Fetch protein sequences from NCBI database
rule fetch_sequences:
    output:
        directory(sequences_dir)
    params:
        bioproject_id=bioproject_id
    shell:
        'python3 fetch_sequences.py --bioproject_id {params.bioproject_id} --output_dir {output}'

# Merge all sequences into one file
rule merge_into_one_file:
    input:
        sequences_dir
    output:
        proteomes_fasta
    shell:
        'python3 merge_into_one_file.py --input_dir {input} --output_file {output}'

# Cluster sequences using MMSeqs2
rule run_clustering:
    input:
        proteomes_fasta
    output:
        directory(mmseqs_results_dir),
        mmseqs_cluster_tsv
    shell:
        'bash clustering.sh {input} {output}'

# Correct clustering to contain orthologs cluster
rule correct_clustering:
    input:
        [sequences_dir, mmseqs_cluster_tsv]
    output:
        directory(clusters_dir)
    shell:
        """
        python3 correct_clustering.py \
        --input_seq_dir {input[0]} \
        --output_dir {output} \
        --cluster_tsv {input[1]}
        """

# Run MSA
rule run_alignment:
    input:
        clusters_dir
    output:
        directory(alignment_dir)
    shell:
        'python3 alignments.py --input_dir {input} --output_dir {output}'

# Infer NJ trees
rule infer_trees:
    input:
        alignment_dir
    output:
        directory(trees_dir)
    shell:
        """
        Rscript infer_trees.R \
            --alignment_dir {input} \
            --output_dir {output}
        """

# Generate consensus
rule consensus_tree:
    input:
        trees_dir
    output:
        directory(consensus_dir)
    shell:
        'python3 consensus_tree.py --trees_dir {input} --output_dir {output}'

# Generate supertree
rule supertree:
    input:
        all_trees_file
    shell:
        'python3 supertree.py --trees_file {input}'

# Infer NJ trees with bootstrap
rule infer_bootstrap_trees:
    input:
        alignment_dir
    output:
        directory(bootstrap_trees_dir)
    params:
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_threshold=bootstrap_threshold
    shell:
        """
        Rscript infer_trees.R \
            --alignment_dir {input} \
            --output_dir {output} \
            --use_bootstrap \
            --bootstrap_replicates {params.bootstrap_replicates} \
            --bootstrap_threshold {params.bootstrap_threshold}
        """

# Generate consensus and supertree for bootstrapped trees
rule consensus_tree_bootstrap:
    input:
        bootstrap_trees_dir
    output:
        directory(bootstrap_consensus_dir)
    shell:
        'python3 consensus_tree.py --trees_dir {input} --output_dir {output}'

rule supertree_bootstrap:
    input:
        bootstrap_trees_file
    shell:
        'python3 supertree.py --trees_file {input}'


# Try another version of clustering correction, keeping paralogs
rule clustering_with_paralogs:
    input:
        [sequences_dir, mmseqs_cluster_tsv]
    output:
        directory(clusters_paralogs_dir)
    shell:
        """
        python3 correct_clustering.py \
        --input_seq_dir {input[0]} \
        --output_dir {output} \
        --cluster_tsv {input[1]}
        """

rule run_alignment_paralogs:
    input:
        clusters_paralogs_dir
    output:
        directory(paralogs_alignment_dir)
    shell:
        'python3 alignments.py --input_dir {input} --output_dir {output}'

rule infer_trees_paralogs:
    input:
        paralogs_alignment_dir
    output:
        directory(paralogs_trees_dir)
    shell:
        """
        Rscript infer_trees.R \
            --alignment_dir {input} \
            --output_dir {output}
        """

rule consensus_tree_paralogs:
    input:
        paralogs_trees_dir
    output:
        directory(paralogs_consensus_dir)
    shell:
        'python3 consensus_tree.py --trees_dir {input} --output_dir {output}'

rule supertree_paralogs:
    input:
        paralogs_trees_file
    shell:
        'python3 supertree.py --trees_file {input}'
