import pandas as pd
import os
from Bio import SeqIO


def get_species_set(input_dir):
    species_set = set()
    for filename in os.listdir(input_dir):
        if filename.endswith('.faa'):
            species = filename[:-4]
            species_set.add(species)
    return species_set

def process_clusters_with_blast(cluster_file, blast_graph_file, output_file, all_species, min_cluster_size=2, keep_paralogs=False):
    # Read processed BLAST data
    blast_results = pd.read_csv(blast_graph_file, sep=' ', header=None)
    blast_results.columns = ['QueryID', 'SubjectID', 'Weight']

    # Read clusters
    clusters = []
    with open(cluster_file, 'r') as f:
        for line in f:
            # Divide line into proteins and extact species names
            proteins = line.strip().split('\t')
            species_to_proteins = {}
            for protein in proteins:
                species = protein.split(':', 1)[0]  # get species name
                # print(species)
                if species not in species_to_proteins:
                    species_to_proteins[species] = []
                species_to_proteins[species].append(protein)
            clusters.append(species_to_proteins)
    print(f'{len(clusters)} clusters in total')

    # Cluster filtering
    filtered_clusters = []
    for cluster in clusters:
        if len(cluster) < min_cluster_size:  # take only large enough clusters
            continue
        
        final_cluster = []
        final_cluster_species = set()

        if keep_paralogs:
            for species, proteins in cluster.items():
                final_cluster_species.add(species)
                final_cluster.extend(proteins)  # add all proteins from given species

        else:  # remove paralogs
            for species, proteins in cluster.items():
                if len(proteins) == 1:
                    final_cluster.append(proteins[0])  # only one protein for species – no paralogs
                    final_cluster_species.add(species)
                else:
                    # More proteins: take the best one from blast results
                    best_protein = None
                    best_score = float('-inf')
                    for protein in proteins:
                        # Sum weights from all comparisons
                        protein_score = blast_results.loc[
                            (blast_results['QueryID'] == protein) | (blast_results['SubjectID'] == protein),
                            'Weight'
                        ].sum()
                        if protein_score > best_score:
                            best_protein = protein
                            best_score = protein_score
                    if best_protein:
                        final_cluster.append(best_protein)
                        final_cluster_species.add(species)
                    else:
                        final_cluster.append(proteins[0])
                        final_cluster_species.add(species)


        # Add cluster to results, if it contains enough proteins
        final_cluster_length = len(final_cluster)
        if final_cluster_length >= min_cluster_size:
            filtered_clusters.append('\t'.join(final_cluster))
            lacking_species = all_species.difference(final_cluster_species)
            print(f'\nAdded cluster {len(filtered_clusters)} of length {final_cluster_length}.\nLacking species: {lacking_species}')

    # Save to file
    with open(output_file, 'w') as f:
        for cluster in filtered_clusters:
            f.write(cluster + '\n')


def clusters_to_fasta(cluster_file, fasta_dir, output_dir):

    os.makedirs(output_dir, exist_ok=True)  # ensure the output directory exists 

    # Read the cluster file
    with open(cluster_file, "r") as file:
        clusters = file.readlines()

    # Process each cluster line
    for cluster_line in clusters:
        proteins = cluster_line.strip().split("\t")

        # Extract the first protein's name (excluding the organism)
        first_protein = proteins[0].split(":")[1]
        output_file = os.path.join(output_dir, f"{first_protein}.fasta")

        # Create a dictionary to hold sequences for this cluster
        cluster_sequences = {}

        for entry in proteins:
            organism, protein = entry.split(":")
            fasta_file = os.path.join(fasta_dir, f"{organism}.fasta")

            if not os.path.exists(fasta_file):
                print(f"Warning: FASTA file for organism {organism} not found.")
                continue

            # Read the corresponding FASTA file and find the protein sequence
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id == protein:
                    # Store the sequence with organism as the header
                    cluster_sequences[organism] = record.seq
                    break

        # Write the cluster sequences to a new FASTA file
        with open(output_file, "w") as output:
            for organism, sequence in cluster_sequences.items():
                output.write(f">{organism}\n{sequence}\n")

    print("FASTA files for clusters have been created in", output_dir)




# species_file = 'input_data/species_ids.txt'
# cluster_file = 'clustering/clusters.txt'
# blast_graph_file = 'blast/blast_graph.mci'

# out_cluster_file = 'clustering/filtered_clusters.txt'
# out_cluster_file_paralogs = 'clustering/filtered_clusters_with_paralogs.txt'
# min_cluster_size = 32

# input_fasta_dir = 'sequences'
# clusters_fasta_dir = 'clusters_fasta'
    
cluster_file = 'clustering/drosophila_clusters.txt'
blast_graph_file = 'blast_drosophila/blast_graph.mci'

out_cluster_file = 'clustering/drosophila_filtered_clusters.txt'
out_cluster_file_paralogs = 'clustering/drosophila_filtered_clusters_with_paralogs.txt'
min_cluster_size = 14

input_fasta_dir = 'drosophila_sequences'
clusters_fasta_dir = 'drosophila_clusters_fasta'

# Lacking species: {'Drosophila_jambulin', 'Drosophila_pandor', 'Drosophila_bipectinat', 'Drosophila_serrat', 'Drosophila_rubid', 'Drosophila_bunnand'}


if __name__ == '__main__':
    # with open(species_file) as sf:
    #     all_species = set(species.strip().replace(' ', '_') for species in sf.readlines())
    all_species = get_species_set(input_fasta_dir)

    # process_clusters_with_blast(cluster_file, blast_graph_file, out_cluster_file, all_species, min_cluster_size)
    process_clusters_with_blast(cluster_file, blast_graph_file, out_cluster_file_paralogs, all_species, min_cluster_size, keep_paralogs=True)
    # clusters_to_fasta(out_cluster_file, input_fasta_dir, clusters_fasta_dir)