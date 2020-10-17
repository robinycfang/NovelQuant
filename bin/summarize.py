#!/usr/bin/env python

import sys
import argparse
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
RI = pd.read_csv(RI, comment = '#', sep = '\t', index_col = 0)
RI = RI.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis = 1)
sum_table = RI.copy()
UJ = pd.read_csv(UJ, comment = '#', sep = '\t')
UJ = UJ.drop(['chromosome', 'start', 'end', 'strand', 'gene'], axis = 1)
UJ_list = set(UJ['transcript'])
for i in UJ_list:
	# if a novel transcript has noth RI and UJ, use RI only
	# if a novel transcript has more than one unique junction, then use the mean of all junction reads
	if i not in RI.index:
		temp = UJ[UJ['transcript'] == i]
		for k in temp.columns:
			if k != 'transcript':
				sum_table.loc[i, k] = temp[k].mean()

# use samtools to extract total sequencing depth of each sample
dep_table = pd.DataFrame(columns = ['depth'])
for line in open(sample_list):
	sample = line.strip('\n')
	print(sample)
	res = subprocess.Popen([st, 'flagstat', sample, '-@', threads], stdout = subprocess.PIPE)
	dep = res.stdout.read()	
	dep = dep.decode('utf-8').split('\n')[0]
	dep = int(dep.split(' ')[0])
	dep_table.loc[sample, 'depth'] = dep
dep_table.to_csv('seq_dep.txt', sep = '\t')

# normalization
for i in sum_table.columns:
	depth = dep_table.loc[i, 'depth']
	sum_table[i] = sum_table[i] * 1000000000 / depth
sum_table.to_csv('NovelQuant_result.txt', sep = '\t')
