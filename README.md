# BLAST presence absence matrix

This script takes the output file from a blast search (format 6, tabular) and converts it into a **presence-absence matrix of 0 and 1**.

The database for the blast search is a multifasta file of reference sequences/genome assemblies (named 'Multifasta.fasta'). The query for the blast search can contain multiple query sequences (example is given of results obtained with a file containing serotype sequences).

The **command used to execute the BLAST search must be the following**: 

```bash
blastn -query Metcalf_serotypes.fasta -db ./Multifasta.fasta -outfmt "6 qseqid sseqid stitle pident qcovs length qlen mismatch gapopen qstart qend sstart send evalue bitscore"-evalue 0.00001 -num_threads 4 -out rel_BLAST_results.txt
```

This is because the script will identify the same columns as in the command above.

## Dependencies

* Biopython v1.76
* pandas v1.0.3

An exaple environment is given as yml file.

## Running the script

Make sure to change the location of your python in the first line of code, after the shebang.
To run the script, simply type:

```bash
presence_absence_matritx.py
```
The script must be run in the same folder where the database (named Multifasta.fasta) and the blast output (named rel_BLAST_results.txt) are stored.

```bash
Optional arguments:
  -h, --help      show this help message and exit
  -pident PIDENT  insert a percentage of identity from 0 to 100 (default: 90)
  -qcov QCOV      insert a percentage of identity from 0 to 100 (default: 80)
 ```
