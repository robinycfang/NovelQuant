#!/usr/bin/env python3

# Find if any complementary introns of annoated genes are retained in the novel transcripts

import pandas as pd
import sys

try:
	complementary_introns = sys.argv[1]
	novel_gtf = sys.argv[2]
except:
	sys.exit('usage: python3 find_retained_introns.py anno_gene_complementary_introns.gtf novel_trans.gtf')

# get complementary introns of annotated genes
anno_gene_union_introns = pd.read_csv(complementary_introns, 
	sep = '\t', chunksize = 10000, header= None, comment = '#',
	names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
anno_gene_union_introns = pd.concat(anno_gene_union_introns, ignore_index = True)

# get all exons of novel transcripts
novel_trans_exons = pd.read_csv(novel_gtf, 
	header = None, sep = '\t', comment = '#',
	names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
novel_trans_exons = novel_trans_exons[novel_trans_exons['type'] == 'exon']

# check if union intronic regions of annotated genes are included in the novel transcripts
retained_introns = ''
# min intronic length to consider
len_filter = 100
for i in anno_gene_union_introns.index:
	in_chr = anno_gene_union_introns.loc[i, 'chromosome']
	in_start = anno_gene_union_introns.loc[i, 'start']
	in_end = anno_gene_union_introns.loc[i, 'end']
	in_strand = anno_gene_union_introns.loc[i, 'strand']

	# select specific chromosome and nearby sites to reduce iteration and speed up the script
	tmp_df = novel_trans_exons[novel_trans_exons['chromosome'] == str(in_chr)]
	tmp_df = tmp_df[tmp_df['start'] - in_start < 10000]
	for k in tmp_df.index:
		ex_strand = tmp_df.loc[k, 'strand']
		# make sure they are both on the same strand
		if in_strand == ex_strand:
			ex_start = tmp_df.loc[k, 'start']
			ex_end = tmp_df.loc[k, 'end']

			# full intron retaintion
			if in_start > ex_start and in_end < ex_end:
				if in_end - in_start >= len_filter:
					retained_introns += str(in_chr) + '\t' + '*' + '\t' + 'intron' + '\t' + str(in_start) + '\t' + str(in_end) + '\t' \
					+ 'IR' + '\t' + in_strand + '\t' + '*' + '\t' + tmp_df.loc[k, 'info'] + '\n'

			# exon 3' extension
			if ex_start < in_start < ex_end and in_end > ex_end:
				if ex_end - in_start >= len_filter:
					retained_introns += str(in_chr) + '\t' + '*' + '\t' + 'intron' + '\t' + str(in_start) + '\t' + str(ex_end) + '\t' \
					+ 'E3E' + '\t' + in_strand + '\t' + '*' + '\t' + tmp_df.loc[k, 'info'] + '\n'

			# exon 5' extension
			if in_start < ex_start and ex_start < in_end < ex_end:
				if in_end - ex_start >= len_filter:
					retained_introns += str(in_chr) + '\t' + '*' + '\t' + 'intron' + '\t' + str(ex_start) + '\t' + str(in_end) + '\t' \
					+ 'E5E' + '\t' + in_strand + '\t' + '*' + '\t' + tmp_df.loc[k, 'info'] + '\n'

			# cassette exon
			if in_start < ex_start and in_end > ex_end:
				if ex_end - ex_start >= len_filter:
					retained_introns += str(in_chr) + '\t' + '*' + '\t' + 'intron' + '\t' + str(ex_start) + '\t' + str(ex_end) + '\t' \
					+ 'CE' + '\t' + in_strand + '\t' + '*' + '\t' + tmp_df.loc[k, 'info'] + '\n'

with open('retained_introns.gtf', 'w') as w:
	w.write(retained_introns)
