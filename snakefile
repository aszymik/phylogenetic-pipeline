# Input/output files and directories

bioproject_id = "PRJNA736147"
sequences_dir = "sequences"
proteomes_fasta = "clustering/proteomes.fasta"
mmseqs_results_dir = "clustering/mmseqs_results"
clusters_dir = "clustering/clusters"
mmseqs_cluster_tsv = f"{mmseqs_results_dir}/mmseqs_cluster.tsv"

rule fetch_sequences:
    output:
        directory=sequences_dir
    params:
        bioproject_id=bioproject_id
    script:
        "fetch_sequences.py"

rule merge_into_one_file:
    input:
        directory=sequences_dir
    output:
        fasta=proteomes_fasta
    script:
        "merge_into_one_file.py"

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
        "correct_clustering.py"


