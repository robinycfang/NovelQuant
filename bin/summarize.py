#!/usr/bin/env python3

import sys
import subprocess
import pandas as pd


try:
	RI = sys.argv[1]
	UJ = sys.argv[2]
	sample_list = sys.argv[3]
	st = sys.argv[4]
	threads = sys.argv[5]
except:
	sys.exit('usage: python summarize.py RI_counts.txt UJ_counts.txt sample_list.txt samtools_path threads')


# merge RI_counts.txt and UJ_counts.txt
# if a novel transcript has both RI and UJ, use RI only
sum_table = pd.DataFrame()
RI_list = []
RI = pd.read_csv(RI, sep = '\t', comment = '#')
for i in RI.index:
	trans = RI.loc[i, 'Geneid']
	RI_list.append(trans)
	for k in RI.columns[6:]:
		sum_table.loc[trans, k] = RI.loc[i, k]

UJ = pd.read_csv(UJ, sep = '\t')
UJ_list = set(UJ['transcript'])
for i in UJ_list:
	if i not in RI_list:
		tmp = UJ[UJ['transcript'] == i].reset_index()
		# for the novel transcripts with only one unique junction
		if len(tmp) == 1:
			for k in tmp.columns[7:]:
				sum_table.loc[i, k] = tmp.loc[0, k]
		# for the novel transcripts with multiple unique junctions, get the means
		else:
			tmp = tmp.iloc[:, 7:].mean().to_frame().T
			for k in tmp.columns:
				sum_table.loc[i, k] = tmp.loc[0, k]
sum_table.to_csv('summary_raw.txt', sep = '\t', index = True)

# use samtools to extract total sequencing depth of each sample
dep_list = {}
dep_holder = ''
for line in open(sample_list):
	sample = line.strip('\n')
	res = subprocess.Popen([st, 'flagstat', sample, '-@', threads], stdout = subprocess.PIPE)
	dep = res.stdout.read()	
	dep = dep.decode('utf-8').split('\n')[0]
	dep = int(dep.split(' ')[0])
	dep_list[sample] = dep
	dep_holder += sample + '\t' + str(dep) + '\n'
with open('sample_depth.txt', 'w') as w:
	w.write(dep_holder)

# normalization
for sample in sum_table.columns:
	dep = dep_list[sample]
	sum_table[sample] = sum_table[sample] * 1000000000 / dep
sum_table.to_csv('summary_norm.txt', sep = '\t', index = True)