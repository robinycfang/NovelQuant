import sys
import argparse
import subprocess
import os

# find complementary introns in each annotated gene
def complement():
	parser = argparse.ArgumentParser(usage = 'python3 NovelQuant complement -i annotation.gtf')
	parser.add_argument('complement')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', help = 'the input gtf file that has information of annotated transcripts', dest = 'i')
	args = parser.parse_args()

	subprocess.call([sys.executable, path + 'complementary_introns.py', args.i])	

# find introns that are retained in the novel transcripts
def findRI():
	parser = argparse.ArgumentParser(usage = 'python3 NovelQuant findRI -i annotation.gtf')
	parser.add_argument('findRI')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', help = 'the gtf file of complementary introns', dest = 'i')
	required.add_argument('-k', help = 'the gtf file of novel transcripts', dest = 'k')
	args = parser.parse_args()

	subprocess.call([sys.executable, path + 'find_retained_introns.py', args.i, args.k])

# quantify novel transcripts by retained introns
def quantRI():
	parser = argparse.ArgumentParser(usage = 'python3 NovelQuant quantRI -i retained_introns.gtf -k BAM -p featureCounts_path')
	parser.add_argument('quantRI')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i', help = 'the gtf file of introns that are retained in novel transcripts', dest = 'i')
	required.add_argument('-k', help = 'The BAM file(s) to process', dest = 'k')
	required.add_argument('-p', help = 'Path to featureCounts if not in the environmental variables', dest = 'p')
	required.add_argument('-t', help = 'Threads to use in featureCounts', dest = 't', type = str)
	args = parser.parse_args()

	subprocess.call([args.p, '-a', args.i, '-o', 'RI_counts.txt', '-t', 'intron', '-g', 'transcript_id', \
			'--minOverlap', '50', '-T', args.t, args.k])


path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
mode = sys.argv[1]

if mode == 'complement':
	complement()

if mode == 'findRI':
	findRI()

if mode == 'quantRI':
	quantRI()
