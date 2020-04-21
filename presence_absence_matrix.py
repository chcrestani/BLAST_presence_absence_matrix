#!/Users/chiara/anaconda3/bin/python

######
# This script takes the output from a blastn output file and converts it into a presence/absence matrix.
# The input file for blast can be a multifasta file (here named "Multifasta.fasta"). The database for blast can also
# contain multiple search sequences.
# The command to execture the BLAST search is:
#
# blastn -query Multifasta.fasta -db ./Capsular_serotype_GBS/Metcalf_serotypes.fasta -outfmt "6 qseqid sseqid stitle 
# pident qcovs length slen mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 0.00001 -num_threads 4 
# -out rel_BLAST_results.txt
#
# This script was written by Chiara Crestani 20200421
#####

# Import libraries
import pandas as pd
from Bio import SeqIO

# Create a dataframe (df) from the tabular blast output file. The header will match the command given for the blast search.

with open('rel_BLAST_results.txt', 'r') as infile: #opening text file for reading (r)
    col_names=['qseqid','sseqid','stitle','pident','qcovs','length','slen','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    df = pd.read_csv(infile, sep='\t', header=None, names=col_names)

    
# Add a column that calculates and gives the query coverage as a percentage

df['qcov'] = df['length'] / df['slen']*100


# Filter the dataframe based on minimum thresholds for percentage of identity and query coverage

df_fil = df[(df['pident'] >= 94) & (df['qcov'] >= 99)] #these will be substituted with variables given with flags


# A dataframe grouping a list of positive results per query sequence id is created.
# A dataframe holding a matrix of presence/absence is crated based on that list.

res_list = df_fil.groupby('qseqid')['sseqid'].apply(list).reset_index()
pa_matrix = res_list.join(pd.get_dummies(res_list['sseqid'].apply(pd.Series).stack()).sum(level=0)).drop('sseqid', 1)


# A list of all query sequence ids in the blast search is extracted from the multifasta file

allids = []
for record in SeqIO.parse("Multifasta.fasta", "fasta"):
    allids.append(record.id)


# A list of query sequence ids from the filtered dataframe of positives is created

resids = df_fil['qseqid'].unique().tolist()


# A list of query sequence ids missing from the filtered dataframe of positives is created and added to the matrix dataframe

s = set(resids)
missing_ids = [x for x in allids if x not in s]

pa_matrix = pa_matrix.append(pd.DataFrame({'qseqid':missing_ids}))


# The matrix dataframe is polished with index reset, filling null vales as 0 and converting floats to integers

pa_matrix.reset_index(drop=True, inplace=True)
pa_matrix.fillna(value=0, inplace=True)

res_col = df_fil['sseqid'].unique().tolist()
pa_matrix[res_col] = pa_matrix[res_col].astype(int)


# The matrix is saved as a csv file 

with open('matrix.csv', 'w') as csv:
    pa_matrix.to_csv(csv, index=False)
