# BLAST database from protein sequences file
makeblastdb -in blast/proteomes.fasta -dbtype prot -out blast/coronavirus_proteomes;

# # outfmt 6: tabular output
# blastp -query blast/proteomes.fasta -db blast/coronavirus_proteomes -outfmt 6 -out blast/blast_results.tsv -evalue 1e-5 -num_threads 4;
blastp -query blast/proteomes.fasta -db blast/coronavirus_proteomes -outfmt 6 -out blast/blast_results.tsv -num_threads 4;

# sudo apt install ncbi-blast+-legacy
# makeblastdb -in blast/proteomes.fasta -dbtype prot -out blast/coronavirus_proteomes;
# blastall -p blastp -i blast/proteomes.fasta -d blast/coronavirus_proteomes -m8 -o blast/blast_results.tsv