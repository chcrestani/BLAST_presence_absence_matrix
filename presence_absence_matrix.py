import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML

with open('rel_BLAST_results.txt', 'r') as infile: #opening text file for reading (r)
    col_names=['qseqid','sseqid','stitle','pident','qcovs','length','slen','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    df = pd.read_csv(infile, sep='\t', header=None, names=col_names)
    
df['qcov'] = df['length'] / df['slen']*100

df_fil = df[(df['pident'] >= 94) & (df['qcov'] >= 99)] #these will be substituted with variables from flag

res_list = df_fil.groupby('qseqid')['sseqid'].apply(list).reset_index()

pa_matrix = res_list.join(pd.get_dummies(res_list['sseqid'].apply(pd.Series).stack()).sum(level=0)).drop('sseqid', 1)


allids = []
for record in SeqIO.parse("Multifasta.fasta", "fasta"):
    allids.append(record.id)

resids = df_fil['qseqid'].unique().tolist()

s = set(resids)
excluded_ids = [x for x in allids if x not in s]

pa_matrix = pa_matrix.append(pd.DataFrame({'qseqid':excluded_ids}))

pa_matrix.reset_index(drop=True, inplace=True)
pa_matrix.fillna(value=0, inplace=True)

res_col = df_fil['sseqid'].unique().tolist()
pa_matrix[res_col] = pa_matrix[res_col].astype(int)


with open('matrix.csv', 'w') as csv:
    pa_matrix.to_csv(csv, index=False)
