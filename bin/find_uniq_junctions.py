# Find unique junctions of novel transcripts

import pandas as pd
import re
import csv
import sys

try:
    eej_gtf = sys.argv[1]  
except:
    sys.exit('usage: find_uniq_junctions.py eej.gtf')

# intron.gtf is too big for pd.read_csv
eej_gtf = pd.read_csv(eej_gtf, sep = '\t', chunksize = 10000, header= None,\
                        names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
eej_gtf = pd.concat(eej_gtf, ignore_index = True)
# convert the dtype of column 'chromosome' to string, otherwise it will have a mix of different dtypes which will cause problems
eej_gtf['chromosome'] = eej_gtf['chromosome'].astype('str')

# append two columns of gene and transcript names
gene = []
trans = []
for entry in eej_gtf['info']:
    each_gene = re.findall(r"gene_id \"(.*?)\";", entry)[0]
    each_trans = re.findall(r"transcript_id \"(.*?)\";", entry)[0]
    gene.append(each_gene)
    trans.append(each_trans)
eej_gtf['gene'] = gene
eej_gtf['transcript'] = trans

# dedup based on chr, start and end to find unique junctions
uniq_junc = eej_gtf.drop_duplicates(subset = ['chromosome', 'start', 'end'], keep = False)
uniq_junc_novel_trans = uniq_junc[-uniq_junc['transcript'].str.contains('ENST')]
uniq_junc_novel_trans = uniq_junc_novel_trans.reset_index(drop = True)
uniq_junc_novel_trans.to_csv('uniq_eej.gtf', index = False, sep = '\t', quoting = csv.QUOTE_NONE)