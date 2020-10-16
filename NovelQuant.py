#!/usr/bin/env python

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
	required.add_argument('-p', help = 'Path to featureCounts if not in the environmental variables', dest = 'p')
	required.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str)
	args = parser.parse_args()

	# quantify novel transcripts by retained introns
	cmd = [args.p, '-a', args.r, '-o', 'RI_counts.txt', '-t', 'intron', '-g', 'transcript_id', '--minOverlap', '50', '-T', args.t]
	for line in open(args.l):
		line = line.strip('\n')
		cmd.append(line)
	subprocess.call(cmd)

def findUJ():
	parser = argparse.ArgumentParser(usage = 'python NovelQuant findUJ -i anno_novel.gtf')
	parser.add_argument('findUJ')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', help = 'The gtf file of exons of both annotated and novel transcripts', dest = 'i')
	args = parser.parse_args()
	
	# extract exon-exon junctions from both annotated and novel transcripts
	subprocess.call([sys.executable, path + 'extract_junctions.py', args.i])
	# find the junctions that are unique to the novel transcripts
	subprocess.call([sys.executable, path + 'find_uniq_junctions.py', 'eej.gtf'])

def quantUJ():
	parser = argparse.ArgumentParser(usage = 'python3 NovelQuant quantUJ -i anno_novel.gtf -e uniq_eej.gtf -l sample_list.txt \
						-p featureCounts_path -t threads')
	parser.add_argument('quantUJ')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', help = 'The gtf file of exons of both annotated and novel transcripts', dest = 'i')	
	required.add_argument('-e', help = 'The gtf file of unique junctions. i.e., the output of findUJ, uniq_eej.gtf', dest = 'e')
	required.add_argument('-l', help = 'A list of BAM file(s) to be processes. \
				Each line should be the path of each bam file.', dest = 'l')
	required.add_argument('-p', help = 'Path to featureCounts if not in the environmental variables', dest = 'p')
	required.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str)
	args = parser.parse_args()

	# use featureCounts to quantify junctions of both annotated and novel transcripts
	cmd = [args.p, '-a', args.i, '-o', 'all_counts', '-J', '-T', args.t]
	for line in open(args.l):
		line = line.strip('\n')
		cmd.append(line)
	subprocess.call(cmd)

	# find reads upon unique junctions of novel transcripts
	all_counts = pd.read_csv('all_counts.jcounts', sep = '\t')
	all_counts = all_counts.drop(['PrimaryGene', 'SecondaryGenes', 'Site1_strand', 'Site2_chr', 'Site2_strand'], axis = 1)
	all_counts = all_counts.rename(columns = {'Site1_chr': 'chromosome', 'Site1_location': 'start', 'Site2_location': 'end'})
	all_counts['start'] = all_counts['start'] + 1
	all_counts['end'] = all_counts['end'] - 1
	uniq_eej = pd.read_csv(args.e, sep = '\t')
	merged = pd.merge(uniq_eej, all_counts, how = 'inner', on = ['chromosome', 'start', 'end'], sort = True)
	merged = merged.drop(['source', 'type', 'none1', 'none2', 'info'], axis = 1)
	merged.to_csv('UJ_counts.txt', index = False, sep = '\t')

path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
mode = sys.argv[1]


if mode == 'findRI':
	findRI()

if mode == 'quantRI':
	quantRI()

if mode == 'findUJ':
	findUJ()

if mode == 'quantUJ':
	quantUJ()
