# BLAST presence absence matrix

This script takes the output file (format 6, tabolar) from a blast search (can be blastn, blastp or blastx) and converts it into a **presence-absence matrix of 0 and 1**.

The database for the blast search is a multifasta file of reference sequences/genome assemblies (here named 'Multifasta.fasta'). The query for the blast search can contain multiple query sequences (example is given of results obtained with a file containing serotype sequences).

The **command used to execute the BLASt search must be the following**: 

'''bash
blastn -query Metcalf_serotypes.fasta -db ./Multifasta.fasta -outfmt "6 qseqid sseqid stitle pident qcovs length qlen mismatch gapopen qstart qend sstart send evalue bitscore"-evalue 0.00001 -num_threads 4 -out rel_BLAST_results.txt'''
