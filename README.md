# BLAST_presence_absence_matrix



# The input file for blast can be a multifasta file (here named
# "Multifasta.fasta"). The database for blast can also
# contain multiple search sequences.
# The command to execture the BLAST search was:
#
# blastn -query Metcalf_serotypes.fasta -db ./Multifasta.fasta
# -outfmt "6 qseqid sseqid stitle pident qcovs length qlen
# mismatch gapopen qstart qend sstart send evalue bitscore"
# -evalue 0.00001 -num_threads 4 -out rel_BLAST_results.txt
