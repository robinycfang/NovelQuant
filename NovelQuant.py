#!/usr/bin/env python3

import sys
import argparse
import subprocess
import os
import pandas as pd


def findRI():
	parser = argparse.ArgumentParser(usage = 'python NovelQuant findRI -a annotated.gtf -n novel.gtf')
	parser.add_argument('findRI')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-a', help = 'The gtf file that has annotated transcript information', dest = 'a')
	required.add_argument('-n', help = 'The gtf file that has novel transcript information', dest = 'n')
	args = parser.parse_args()

	# find complementary introns in each annotated gene
	subprocess.call([sys.executable, path + 'complementary_introns.py', args.a])
	# find introns that are retained in the novel transcripts
	subprocess.call([sys.executable, path + 'find_retained_introns.py', 'anno_gene_complementary_introns.gtf', args.n])

def quantRI():
	parser = argparse.ArgumentParser(usage = 'python NovelQuant quantRI -r retained_introns.gtf \
						-l sample_list.txt -p featureCounts_path -t threads')
	parser.add_argument('quantRI')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-r', help = 'The gtf file of introns that are retained in novel transcripts.\
				i.e., the output of findRI, retained_introns.gtf', dest = 'r')
	required.add_argument('-l', help = 'A list of BAM file(s) to be processed. \
				Each line should be the path of each bam file.', dest = 'l')
	required.add_argument('-p', help = 'Path to featureCounts if not in the environmental variables', dest = 'p', default = 'featureCounts')
	required.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str, default = '1')
	args = parser.parse_args()

	# quantify novel transcripts by retained introns
	cmd = [args.p, '-a', args.r, '-o', 'RI_counts.txt', '-t', 'intron', '-g', 'transcript_id', '--minOverlap', '50', '-T', args.t]
	for line in open(args.l):
		line = line.strip('\n')
		cmd.append(line)
	subprocess.call(cmd)

def findUJ():
	parser = argparse.ArgumentParser(usage = 'python NovelQuant findUJ -a annotated.gtf -n novel.gtf')
	parser.add_argument('findUJ')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-a', help = 'The gtf file of exons of annotated transcripts', dest = 'a')
	required.add_argument('-n', help = 'The gtf file of exons of novel transcripts', dest = 'n')
	args = parser.parse_args()

	# extract exon-exon junctions from both annotated and novel transcripts
	subprocess.call([sys.executable, path + 'extract_junctions.py', args.a, args.n])
	# find the junctions that are unique to the novel transcripts
	subprocess.call([sys.executable, path + 'find_uniq_junctions.py', 'eej.gtf'])

def quantUJ():
	parser = argparse.ArgumentParser(usage = 'python3 NovelQuant quantUJ -a annotated.gtf -n novel.gtf -e uniq_eej.gtf -l sample_list.txt \
						-p featureCounts_path -t threads')
	parser.add_argument('quantUJ')
	required = parser.add_argument_group('required arguments')	
	# required.add_argument('-a', help = 'The gtf file of exons of annotated transcripts', dest = 'a')
	required.add_argument('-n', help = 'The gtf file of exons of novel transcripts', dest = 'n')
	required.add_argument('-e', help = 'The gtf file of unique junctions. i.e., the output of findUJ, uniq_eej.gtf', dest = 'e')
	required.add_argument('-l', help = 'A list of BAM file(s) to be processed. \
				Each line should be the path of each bam file.', dest = 'l')
	required.add_argument('-p', help = 'Path to featureCounts if not in the environmental variables', dest = 'p', default = 'featureCounts')
	required.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str, default = '1')
	args = parser.parse_args()

	sample_num = subprocess.Popen(['wc', '-l', args.l], stdout = subprocess.PIPE)
	sample_num = int(sample_num.stdout.read().decode('utf-8'))
	if sample_num <= 10:
		cmd = [args.p, '-a', args.n, '-o', 'novel_counts', '-J', '-T', args.t]
		for line in open(args.l):
			line = line.strip('\n')
			cmd.append(line)
		subprocess.call(cmd)
	else:
		cmd = [args.p, '-a', args.n, '-o', 'novel_counts_1', '-J', '-T', args.t]
		counter = 0
		surfix_name = 1
		for line in open(args.l):
			line = line.strip('\n')
			cmd.append(line)
			counter += 1
			if counter % 10 == 0 or counter == sample_num:
				subprocess.call(cmd)
				surfix_name += 1
				cmd = [args.p, '-a', args.n, '-o', 'novel_counts_' + str(surfix_name), '-J', '-T', args.t]

	uniq_eej = pd.read_csv(args.e, sep = '\t')
	uniq_eej = uniq_eej.drop(['source', 'type', 'none1', 'none2', 'info'], axis = 1)
	merged = uniq_eej.copy()
	for i in os.listdir():
		if 'jcounts' in i:
			all_counts = pd.read_csv(i, sep = '\t')
			all_counts = all_counts.drop(['PrimaryGene', 'SecondaryGenes', 'Site1_strand', 'Site2_chr', 'Site2_strand'], axis = 1)
			all_counts = all_counts.rename(columns = {'Site1_chr': 'chromosome', 'Site1_location': 'start', 'Site2_location': 'end'})
			all_counts['start'] = all_counts['start'] + 1
			all_counts['end'] = all_counts['end'] - 1
			merged = pd.merge(merged, all_counts, how = 'left', on = ['chromosome', 'start', 'end'], sort = True)
	merged = merged.fillna(0)
	merged.to_csv('UJ_counts.txt', index = False, sep = '\t')
	# remove redundant featureCounts outputs (using shell = True to avoid the error caused by qutations)
	subprocess.call('rm novel_counts_*', shell = True)

def summarize():
	parser = argparse.ArgumentParser(usage = 'python NovelQuant sum')
	parser.add_argument('sum')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-r', help = 'Output from quantRI, RI_counts.txt.', dest = 'r')
	required.add_argument('-u', help = 'Output from quantUJ, UJ_counts.txt.', dest = 'u')
	required.add_argument('-l', help = 'A list of BAM file(s) to be processed. \
				Each line should be the path of each bam file.', dest = 'l')
	required.add_argument('-p', help = 'Path to samtools if not in the environmental variables', dest = 'p', default = 'samtools')
	required.add_argument('-t', help = 'Threads to use in samtools.', dest = 't', type = str, default = '1')	
	args = parser.parse_args()

	# merge RI_counts.txt and UJ_counts.txt
	# if a novel transcript has both RI and UJ, use RI only
	sum_table = pd.DataFrame()
	RI_list = []
	RI = pd.read_csv(args.r, sep = '\t', comment = '#')
	for i in RI.index:
		trans = RI.loc[i, 'Geneid']
		RI_list.append(trans)
		for k in RI.columns[6:]:
			sum_table.loc[trans, k] = RI.loc[i, k]

	UJ = pd.read_csv(args.u, sep = '\t')
	UJ_list = set(UJ['transcript'])
	for i in UJ_list:
		if i not in RI_list:
			tmp = UJ[UJ['transcript'] == i].reset_index()
			# for the novel transcripts with only one unique junction
			if len(tmp) == 1:
				for k in tmp.columns[6:]:
					sum_table.loc[i, k] = tmp.loc[0, k]
			# for the novel transcripts with multiple unique junctions, get the means
			else:
				tmp = tmp.iloc[:, 6:].mean().to_frame().T
				for k in tmp.columns:
					sum_table.loc[i, k] = tmp.loc[0, k]
	sum_table.to_csv('summary_raw.txt', sep = '\t', index = True)

	# use samtools to extract total sequencing depth of each sample
	dep_list = {}
	dep_holder = ''
	for line in open(args.l):
		sample = line.strip('\n')
		res = subprocess.Popen(['samtools', 'flagstat', sample, '-@', args.t], stdout = subprocess.PIPE)
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


path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'

if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
	sys.exit('usage: python NovelQuant.py <mode>' + '\n' + \
			'modes: findRI, quantRI, findUJ, quantUJ, sum')
else:
	mode = sys.argv[1]

if mode == 'findRI':
	findRI()
if mode == 'quantRI':
	quantRI()
if mode == 'findUJ':
	findUJ()
if mode == 'quantUJ':
	quantUJ()
if mode == 'sum':
	summarize()