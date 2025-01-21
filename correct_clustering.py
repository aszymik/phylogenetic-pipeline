import os
import subprocess
from Bio import SeqIO
from collections import defaultdict
from functools import lru_cache


def parse_header(header):
    """
    Parses a FASTA header in the format '>species:protein'.

    Args:
        header (str): FASTA header.

    Returns:
        tuple: Species and protein ID.
    """
    species, protein = header.split(":")
    return species.strip(), protein.strip()


def parse_mmseqs_outputs(file_path):
    """
    Parses MMseqs2 clustering output to create a nested dictionary structure.

    Args:
        file_path (str): Path to the MMseqs2 clustering output file.

    Returns:
        dict: Dictionary {cluster_id: {species: list of proteins}}.
    """
    cluster_helper = {}
    with open(file_path, 'r') as file:
        for line in file:
            leader, seq_id = line.strip().split()  # First sequence in cluster is cluster ID
            cluster_helper[leader] = cluster_helper.get(leader, [])
            cluster_helper[leader].append(seq_id)

    clustering = {}
    for cluster_header, seq_header_list in cluster_helper.items():
        _, cluster_id = parse_header(cluster_header)
        this_cluster = {}

        for seq_header in seq_header_list:
            species, protein = parse_header(seq_header)
            this_cluster[species] = this_cluster.get(species, [])
            this_cluster[species].append(protein)

        clustering[cluster_id] = this_cluster

    return clustering


@lru_cache(maxsize=None)
def get_protein_sequence(input_seq_dir, species, protein_tuple):
    """
    Extracts protein sequences from a given FASTA file for specified protein IDs.

    Args:
        input_seq_dir (str): Path to the directory containing all input FASTA files.
        species (str): Species name the proteins come from.
        protein_tuple (tuple): Tuple of protein IDs to extract.

    Returns:
        list: List of protein sequences corresponding to the given IDs.
    """
    sequences = {}
    input_fasta = f'{input_seq_dir}/{species}.faa'
    with open(input_fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_id = record.id.split()[0]  # Extract the first token before the whitespace
            if protein_id in protein_tuple:
                sequences[protein_id] = str(record.seq)

    return [sequences[protein_id] for protein_id in protein_tuple if protein_id in sequences]


def get_best_protein(input_seq_dir, species, proteins_to_choose, rest_of_cluster):
    """
    Performs BLASTP to determine the best protein among candidates based on cluster context.

    Args:
        input_seq_dir (str): Path to the directory containing all input FASTA files.
        proteins_to_choose (list): List of candidate protein IDs.
        rest_of_cluster (dict): Dict of all other species in the cluster with their protein IDs.

    Returns:
        str: The protein ID with the best BLASTP score.
    """
    tmp_dir = "blast_tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    # Prepare input files for BLAST
    proteins_seq = get_protein_sequence(input_seq_dir, species, tuple(proteins_to_choose))
    proteins_fasta = os.path.join(tmp_dir, "proteins_to_choose.fasta")
    with open(proteins_fasta, 'w') as f:
        for protein_id, seq in zip(proteins_to_choose, proteins_seq):
            f.write(f">{protein_id}\n{seq}\n")

    rest_cluster_seq = []
    for other_species, protein_list in rest_of_cluster.items():
        rest_cluster_seq.extend(get_protein_sequence(input_seq_dir, other_species, tuple(protein_list)))
    rest_fasta = os.path.join(tmp_dir, "rest_of_cluster.fasta")
    with open(rest_fasta, 'w') as f:
        for i, seq in enumerate(rest_cluster_seq):
            f.write(f">seq{i}\n{seq}\n")

    # Perform BLASTP
    blast_output = os.path.join(tmp_dir, "blast_results.txt")
    subprocess.run([
        "blastp",
        "-query", proteins_fasta,
        "-subject", rest_fasta,
        "-out", blast_output,
        "-outfmt", "6 qseqid bitscore"
    ])

    # Parse BLAST results
    scores = defaultdict(float)
    with open(blast_output, 'r') as f:
        for line in f:
            query_id, bitscore = line.strip().split("\t")
            scores[query_id] += float(bitscore)

    # Return the protein with the highest total score
    best_protein = max(scores, key=scores.get)

    # Cleanup temporary files
    for tmp_file in [proteins_fasta, rest_fasta, blast_output]:
        os.remove(tmp_file)

    return best_protein


def filter_by_cluster_size(clustering, min_cluster_size):
    """
    Filters clusters based on their size.

    Args:
        clustering (dict): Dictionary of clusters.
        min_cluster_size (int): Minimum number of species required in a cluster.

    Returns:
        dict: Filtered clustering dictionary.
    """
    return {cluster_id: cluster for cluster_id, cluster in clustering.items() if len(cluster) >= min_cluster_size}


def process_clusters(cluster_tsv, input_seq_dir, output_dir, min_cluster_size=15, remove_paralogs=True):
    """
    Processes clusters, optionally removing paralogs and saving each cluster to a separate file.

    Args:
        cluster_tsv (str): Path to the MMseqs2 clustering output TSV file.
        input_fasta (str): Path to the input FASTA file with all sequences.
        output_dir (str): Directory to save processed cluster files.
        min_cluster_size (int): Minimum size of clusters to retain.
        remove_paralogs (bool): Whether to remove paralogs within clusters.
    """
    os.makedirs(output_dir, exist_ok=True)  # Create output directory if not exists

    clusters = parse_mmseqs_outputs(cluster_tsv)
    clusters = filter_by_cluster_size(clusters, min_cluster_size)

    for cluster_id, cluster in clusters.items():
        if remove_paralogs:
            for species, protein_list in cluster.items():
                if len(protein_list) > 1:
                    rest_of_cluster = {}
                    for other_species, other_proteins in cluster.items():
                        if other_species != species:
                            rest_of_cluster[other_species] = other_proteins

                    chosen_protein = get_best_protein(input_seq_dir, species, protein_list, rest_of_cluster)
                    cluster[species] = [chosen_protein]

        # Save the cluster to a separate file
        output_file = os.path.join(output_dir, f"{cluster_id}.fasta")

        with open(output_file, 'w') as out_fasta:
            for species, protein_list in cluster.items():
                sequences = get_protein_sequence(input_seq_dir, species, tuple(protein_list))
                for seq_id, seq in zip(protein_list, sequences):
                    out_fasta.write(f">{seq_id}\n{seq}\n")


if __name__ == '__main__':
    input_seq_dir = 'sequences'
    input_fasta = "clustering/proteomes.fasta"
    mmseq_output_dir = "clustering/mmseqs_results"
    cluster_dir = "clustering/clusters"
    cluster_tsv = f'{mmseq_output_dir}/mmseqs_cluster.tsv'

    process_clusters(cluster_tsv, input_seq_dir, cluster_dir, remove_paralogs=True)

