# Phylogenetic pipeline

This pipeline automates fetching, merging, clustering, correcting protein sequences, and inferring trees for phylogenetic analysis.

## Overview

The pipeline follows these main steps:

1. **Fetch sequences**: Download protein FASTA files from a given BioProject ID (NCBI).
2. **Merge sequences**: Combine all downloaded `.faa` files into a single FASTA file, modifying headers as needed.
3. **Clustering**: Run MMseqs2 (easy-linclust) on the merged sequences to cluster them based on sequence similarity.
4. **Correct clustering**: Filter clusters, remove paralogs, and generate final cluster files suitable for downstream phylogenetic analysis.
5. **Multi-sequence alignment**: Sequences within clusters are aligned using MAFFT.
6. **Gene tree construction**: Neighbor-joining trees are constructed for each family using Phangorn.
7. **Consensus tree generation**: Consensus trees are generated using IQ-TREE.
8. **Supertree generation**: Supertrees are generated using Fasturec.

---

## System requirements

Before starting, make sure the following tools are installed:
- [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/mafft_7.xx-with-extensions-src.tgz)
- [Fasturec](https://bitbucket.org/pgor17/fasturec.git)


### Instaling NCBI BLAST+
```bash
# Linux (Debian/Ubuntu)
sudo apt install ncbi-blast+

# or manually:
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*-x64-linux.tar.gz
tar -xzf ncbi-blast-*-x64-linux.tar.gz
export PATH=$PATH:/path/to/blast/bin
```

### Installing MAFFT
```bash
# Linux (Debian/Ubuntu)
sudo apt install mafft

# or manually:
wget https://mafft.cbrc.jp/alignment/software/mafft_7.xx-with-extensions-src.tgz
tar zxvf mafft_7.xx-with-extensions-src.tgz
cd mafft-7.xx-with-extensions/core/
sudo make install
```

### Installing Fasturec

```bash
git clone https://bitbucket.org/pgor17/fasturec.git
cd fasturec || exit
make
```

All Python packages that have to be installed with `conda` are listed in `conda_requirements.txt` and can be installed using `conda create -n env --file package-list.txt`

All other Python packages are listed in `requirements.txt` and can be installed using `pip install -r requirements.txt`.


## Running the pipeline

When all dependencies are installed, you can run the pipeline using command

```bash
snakemake --cores 1
```

Adjust `--cores` based on your available computing resources.

By default, the `snakefile` uses:
```
bioproject_id = "PRJNA736147"
sequences_dir = "sequences"
proteomes_fasta = "clustering/proteomes.fasta"
mmseqs_results_dir = "clustering/mmseqs_results"
clusters_dir = "clustering/clusters"
mmseqs_cluster_tsv = "clustering/mmseqs_results/mmseqs_cluster.tsv"
```
You can modify these variables in the `snakefile` to match your data paths.


## Pipeline Steps & Scripts

Below is a summary of each script and its function in the workflow – although they form a pipeline, each of them can be used as a standalone script.

1. `fetch_sequences.py`

Fetches protein sequences from a given BioProject (NCBI) using Biopython’s Entrez module.
Downloads .faa files (protein FASTA) for each assembly in the BioProject.

Usage Example:

```
python fetch_sequences.py \
    --bioproject_id PRJNA736147 \
    --output_dir sequences
```

2. `merge_into_one_file.py`

Merges all .faa files in a directory into a single FASTA file (by default, `proteomes.fasta`).
Reformats headers to include species name.

Usage Example:

```
python merge_into_one_file.py \
    --input_dir sequences \
    --output_file clustering/proteomes.fasta
```

3. `clustering.sh`

Runs MMseqs2 (easy-linclust) on the merged FASTA file.
Outputs the clustering results in a specified directory.

Usage Example:

```
bash clustering.sh clustering/proteomes.fasta clustering/mmseqs_results
```

4. `correct_clustering.py`

Parses MMseqs2 clustering TSV file.
Removes paralogs by selecting the best protein via BLASTP-based scoring against the rest of the cluster.
Filters clusters by minimum cluster size.
Writes each cluster’s final set of protein sequences into separate FASTA files.

Usage Example:

```
python correct_clustering.py \
    --input_seq_dir sequences \
    --output_dir clustering/clusters \
    --cluster_tsv clustering/mmseqs_results/mmseqs_cluster.tsv \
    --min_cluster_size 15 \
    --remove_paralogs True
```

5. `alignments.py`

Aligns sequences within clusters using MAFFT. Outputs aligned FASTA files for each cluster.

Usage Example:
```
python alignments.py \
    --input_dir clustering/clusters \
    --output_dir msa
```

6. `infer_trees.R`

Infers phylogenetic trees using the alignments as input. Utilizes Phangorn to construct neighbor-joining trees.

Usage Example:
```
Rscript infer_trees.R \
    --alignment_dir msa \
    --output_dir trees
```

7. `consensus_tree.py`

Generates greedy and majority consensus trees for the inferred trees using IQ-TREE.

Usage Example:
```
python consensus_tree.py \
    --trees_dir trees \
    --output_dir consensus
```

8. `supertree.py`
Generates supertrees from the consensus trees using Fasturec. Outputs a single supertree for all families.

Usage Example:
```
python supertree.py \
    --trees_file consensus/all_trees.treefile
```

