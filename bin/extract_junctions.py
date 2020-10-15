# Find the exon-exon junctions of transcripts in a gtf file.

import pandas as pd
import re
import sys

try:
	gtf_path = sys.argv[1]	# this gtf file must contain both annotated and novel transcripts
except:
	sys.exit('usage: extract_junctions.py anno_novel_trans.gtf')

gtf = pd.read_csv(gtf_path, sep = '\t', header = None, comment = '#')
exons = gtf[gtf[2] == 'exon']
# create the last column to indicate transcript names
transcripts = []
for i in exons[8]:
    trans_name = re.findall(r"transcript_id \"(.*?)\";", i)[0]
    transcripts.append(trans_name)
exons[9] = transcripts
# sort the exon table by chromosome, transcript name and then the start position
exons = exons.sort_values(by = [0, 9, 3]).reset_index(drop = True)

# iterate the gtf file to find exon-exon junctions 
mul = open('eej.gtf', 'a')
meter = 1

try:
    for i in range(len(exons)):
        chromosome = exons.loc[i, 0]
        source = exons.loc[i, 1]
        start_exon = exons.loc[i, 3]
        end_exon = exons.loc[i, 4]
        strand = exons.loc[i, 6]
        info = exons.loc[i, 8]
        this_trans = exons.loc[i, 9]
        next_trans = exons.loc[i + 1, 9]

        # for multiple-exon transcripts
        if this_trans == next_trans:
            start_eej = end_exon + 1
            end_eej = exons.loc[i + 1, 3] - 1
            multi_ex_trans = str(chromosome) + '\t' + source + '\t' + 'eej' + '\t' + str(start_eej) + '\t' + str(end_eej) + \
                                '\t' + '.' + '\t' + strand + '\t' + '.' + '\t' + info + '\n'
            mul.write(multi_ex_trans)

        print(str(meter) + '/' + str(len(exons)) + ' done.')
        meter += 1

# last line causes an error
except KeyError:
        print('All done.')

mul.close()