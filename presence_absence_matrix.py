#!/Users/chiara/anaconda3/bin/python

######
# This script takes the output from a blastn output file and converts
#  it into a presence/absence matrix.
# The input file for blast can be a multifasta file (here named
# "Multifasta.fasta"). The database for blast can also
# contain multiple search sequences.
# The command to execture the BLAST search is:
#
# blastn -query Multifasta.fasta -db ./Metcalf_serotypes.fasta
# -outfmt "6 qseqid sseqid stitle pident qcovs length slen mismatch
# gapopen qstart qend sstart send evalue bitscore" -evalue 0.00001
# -num_threads 4 -out rel_BLAST_results.txt
######

"""
@author: Chiara Crestani (chcrestani)
"""


# Import libraries
import pandas as pd
from Bio import SeqIO
import argparse

# Script flags and messages
parser = argparse.ArgumentParser(description='Convert the tabular \
                                 output (outfmt 6) from a blastn \
                                 search into a presence/absence matrix.')
parser.add_argument('-pident', default=90, type=int,
                    help='insert a percentage of identity \
                    from 0 to 100 (default: 90)')
parser.add_argument('-qcov', default=80, type=int,
                    help='insert a percentage of identity \
                    from 0 to 100 (default: 80)')

args = parser.parse_args()

pident = args.pident
qcov = args.qcov

# Create a dataframe (df) from the tabular blast output file.
# The header will match the command given for the blast search.
col_names=['qseqid','sseqid','stitle','pident',
           'qcovs','length','slen','mismatch',
           'gapopen','qstart','qend','sstart',
           'send','evalue','bitscore']
df = pd.read_csv('rel_BLAST_results.txt',
                 sep='\t', header=None,
                 names=col_names)

# Add a column that calculates and gives the query coverage
# as a percentage
df['qcov'] = df['length'] / df['slen']*100

# Filter the dataframe based on minimum thresholds for percentage
# of identity and query coverage
df_fil = df[(df['pident'] >= pident) & (df['qcov'] >= qcov)]

# A dataframe grouping a list of positive results per query
# sequence id is created.
res_list = df_fil.groupby('qseqid')['sseqid']\
                          .apply(list).reset_index()

# A dataframe holding a matrix of presence/absence is created
# based on that list.
pa_matrix = res_list.join(pd.get_dummies(res_list['sseqid']\
                                         .apply(pd.Series)\
                                         .stack())\
                                         .sum(level=0))\
                          .drop('sseqid', 1)

# A list of all query sequence ids in the blast search is
# extracted from the multifasta file
allids = []
for record in SeqIO.parse("Multifasta.fasta", "fasta"):
    allids.append(record.id)

# A list of query sequence ids from the filtered dataframe of
# positives is created
resids = df_fil['qseqid'].unique().tolist()

# A list of query sequence ids missing from the filtered
# dataframe of positives is created and added to the matrix
# dataframe
s = set(resids)
missing_ids = [x for x in allids if x not in s]

pa_matrix = pa_matrix.append(pd.DataFrame({'qseqid':missing_ids}))

# The matrix dataframe is polished with index reset,
# filling null vales as 0 and converting floats to integers
pa_matrix.reset_index(drop=True, inplace=True)
pa_matrix.fillna(value=0, inplace=True)

res_col = df_fil['sseqid'].unique().tolist()
pa_matrix[res_col] = pa_matrix[res_col].astype(int)


# The matrix is saved as a csv file
pa_matrix.to_csv('matrix.csv', index=False)
