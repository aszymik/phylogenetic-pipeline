# BLAST database from file with protein sequences
makeblastdb -in proteomes.fasta -dbtype prot -out coronavirus_proteomes;

# outfmt 6 – tabular output
blastp -query proteomes.fasta -db coronavirus_proteomes -outfmt 6 -out blast_results.tsv -evalue 1e-5 -num_threads 4;

