import pandas as pd


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


species_file = 'input_data/species_ids.txt'
cluster_file = 'clustering/clusters.txt'
blast_graph_file = 'blast/blast_graph.mci'
output_file = 'clustering/filtered_clusters.txt'
output_file_paralogs = 'clustering/filtered_clusters_with_paralogs.txt'
min_cluster_size = 32

if __name__ == '__main__':
    with open(species_file) as sf:
        all_species = set(species.strip().replace(' ', '_') for species in sf.readlines())

    process_clusters_with_blast(cluster_file, blast_graph_file, output_file, all_species, min_cluster_size)
    process_clusters_with_blast(cluster_file, blast_graph_file, output_file_paralogs, all_species, min_cluster_size, keep_paralogs=True)
