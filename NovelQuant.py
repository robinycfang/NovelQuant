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
	parser = argparse.ArgumentParser(usage = 'python NovelQuant quantRI -r retained_introns.gtf -l sample_list.txt')
	parser.add_argument('quantRI')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-r', help = 'The gtf file of introns that are retained in novel transcripts.\
				i.e., the output of findRI, retained_introns.gtf', dest = 'r')
	required.add_argument('-l', help = 'A list of BAM file(s) to be processed. \
				Each line should be the path of each bam file.', dest = 'l')
	parser.add_argument('-p', help = 'Path to featureCounts if not in the environmental variables', dest = 'p', default = 'featureCounts')
	parser.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str, default = '1')
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
	parser = argparse.ArgumentParser(usage = 'python3 NovelQuant quantUJ -n novel.gtf -e uniq_eej.gtf -l sample_list.txt')
	parser.add_argument('quantUJ')
	required = parser.add_argument_group('required arguments')	
	required.add_argument('-n', help = 'The gtf file of exons of novel transcripts', dest = 'n')
	required.add_argument('-e', help = 'The gtf file of unique junctions. i.e., the output of findUJ, uniq_eej.gtf', dest = 'e')
	required.add_argument('-l', help = 'A list of BAM file(s) to be processed. \
				Each line should be the path of each bam file.', dest = 'l')
	parser.add_argument('-p', help = 'Path to featureCounts if not in the environmental variables', dest = 'p', default = 'featureCounts')
	parser.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str, default = '1')
	args = parser.parse_args()

	sample_num = subprocess.Popen(['wc', '-l', args.l], stdout = subprocess.PIPE)
	sample_num = int(sample_num.stdout.read().decode('utf-8').split(' ')[0])
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
	subprocess.call('rm novel_counts*', shell = True)

def summarize():
	parser = argparse.ArgumentParser(usage = 'python NovelQuant sum -r RI_counts.txt -u UJ_counts.txt -l sample_list.txt')
	parser.add_argument('sum')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-r', help = 'Output from quantRI, RI_counts.txt.', dest = 'r')
	required.add_argument('-u', help = 'Output from quantUJ, UJ_counts.txt.', dest = 'u')
	required.add_argument('-l', help = 'A list of BAM file(s) to be processed. \
				Each line should be the path of each bam file.', dest = 'l')
	parser.add_argument('-st', help = 'Path to samtools if not in the environment variables', dest = 'st', default = 'samtools')
	parser.add_argument('-a', help = 'The gtf file that has annotated transcript information', dest = 'a')
	parser.add_argument('-n', help = 'The gtf file that has novel transcript information', dest = 'n')
	parser.add_argument('-fc', help = 'Path to featureCounts if not in the environmental variables', dest = 'fc', default = 'featureCounts')
	parser.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str, default = '1')

	args = parser.parse_args()

	subprocess.call([sys.executable, path + 'summarize.py', args.r, args.u, args.l, args.st, args.a, args.n, args.fc, args.t])


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