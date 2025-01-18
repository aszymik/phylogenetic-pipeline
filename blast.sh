# BLAST database from protein sequences file
makeblastdb -in blast/proteomes.fasta -dbtype prot -out blast/coronavirus_proteomes;

# outfmt 6: tabular output
blastp -query blast/proteomes.fasta -db blast/coronavirus_proteomes -outfmt 6 -out blast/blast_results.tsv -evalue 1e-5 -num_threads 4;

