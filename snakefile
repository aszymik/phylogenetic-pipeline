# Input/output files and directories

bioproject_id = 'PRJNA736147'
sequences_dir = 'sequences'
proteomes_fasta = 'clustering/proteomes.fasta'
mmseqs_results_dir = 'clustering/mmseqs_results'
clusters_dir = 'clustering/clusters'
mmseqs_cluster_tsv = f'{mmseqs_results_dir}/mmseqs_cluster.tsv'
alignment_dir = 'msa'
trees_dir = 'trees'
consensus_dir = 'consensus'

rule fetch_sequences:
    output:
        directory=sequences_dir
    params:
        bioproject_id=bioproject_id
    script:
        'fetch_sequences.py'

rule merge_into_one_file:
    input:
        directory=sequences_dir
    output:
        fasta=proteomes_fasta
    script:
        'merge_into_one_file.py'

rule run_clustering:
    input:
        fasta=proteomes_fasta
    output:
        directory=mmseqs_results_dir
    shell:
        """
        bash clustering.sh {input.fasta} {output.directory}
        """

rule correct_clustering:
    input:
        seq_dir=sequences_dir,
        cluster_tsv=mmseqs_cluster_tsv
    output:
        directory=clusters_dir
    script:
        'correct_clustering.py'

rule run_alignment:
    input:
        directory=clusters_dir
    output:
        directory=alignment_dir
    script:
        'alignments.py'

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

