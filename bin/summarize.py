#!/usr/bin/env python3

import sys
import subprocess
import pandas as pd
import re
from statistics import median

pd.options.mode.chained_assignment = None

######################################### merge RI_counts.txt and UJ_counts.txt

# get parameters for merging
try:
    RI = sys.argv[1]
    UJ = sys.argv[2]
    sample_list = sys.argv[3]
    st = sys.argv[4]
except:
    sys.exit('usage: python summarize.py RI_counts.txt UJ_counts.txt sample_list.txt \
            (provide anno.gtf and novel.gtf for calculating expression percentages of novel transcripts)')

# if a novel transcript has both RI and UJ, use RI only
sum_table = pd.DataFrame(columns = ['quant_type'])
RI_list = []
RI = pd.read_csv(RI, sep = '\t', comment = '#', index_col = 0)
for trans in RI.index:
    RI_list.append(trans)
    for k in RI.columns[5:]:
        sum_table.at[trans, k] = RI.at[trans, k]
sum_table['quant_type'] = 'retained_intron'

UJ = pd.read_csv(UJ, sep = '\t')
UJ_list = set(UJ['transcript'])
for i in UJ_list:
    if i not in RI_list:
        tmp = UJ[UJ['transcript'] == i].reset_index()
        # for the novel transcripts with only one unique junction
        if len(tmp) == 1:
            for k in tmp.columns[7:]:
                sum_table.at[i, k] = tmp.at[0, k]
        # for the novel transcripts with multiple unique junctions, get the means
        else:
            tmp = tmp.iloc[:, 7:].mean().to_frame().T
            for k in tmp.columns:
                sum_table.at[i, k] = tmp.at[0, k]
sum_table['quant_type'].fillna('unique_junction', inplace = True)
#sum_table.to_csv('summary_raw.txt', sep = '\t', index = True)

# use samtools to extract total sequencing depth of each sample
dep_list = {}
dep_holder = ''
for line in open(sample_list):
    sample = line.strip('\n')
    res = subprocess.Popen([st, 'flagstat', sample], stdout = subprocess.PIPE)
    dep = res.stdout.read() 
    dep = dep.decode('utf-8').split('\n')[0]
    dep = int(dep.split(' ')[0])
    dep_list[sample] = dep
    dep_holder += sample + '\t' + str(dep) + '\n'
with open('sample_depth.txt', 'w') as w:
    w.write(dep_holder)

# normalization
for sample in sum_table.columns[1:]:
    dep = dep_list[sample]
    sum_table[sample] = sum_table[sample] * 1000000000 / dep
sum_table.to_csv('NovelQuant_final.txt', sep = '\t', index = True, index_label = 'transcript')

########################################### calculate expression percentage

if len(sys.argv) == 9:
    anno_gtf = sys.argv[5]
    novel_gtf = sys.argv[6]
    fc = sys.argv[7]
    threads = sys.argv[8]

    # a function to get transcript and gene name matches from GTF files
    def find_match(in_gtf):
        match_df = pd.DataFrame(columns = ['transcript', 'gene'])
        in_gtf = pd.read_csv(in_gtf, comment = '#', sep = '\t', header = None, low_memory = False)
        # get only the transcript entries; watch out flair gtf output has no transcript entries
        in_gtf = in_gtf[in_gtf.iloc[:, 2] == 'transcript']
        # working on the series so doesn't need to specific axis for apply func
        match_df['transcript'] = in_gtf.iloc[:, -1].apply(lambda x: re.findall(r'transcript_id "(.*?)"', x)[0])
        match_df['gene'] = in_gtf.iloc[:, -1].apply(lambda x: re.findall(r'gene_id "(.*?)"', x)[0])
        match_df.set_index('transcript', inplace = True)

        return match_df

    # use it
    anno_match = find_match(anno_gtf)
    novel_match = find_match(novel_gtf)
    # concate to make the final table
    all_match = pd.concat([anno_match, novel_match], axis = 0)

    # use featureCounts to count annotated transcripts
    cmd = [fc, '-a', anno_gtf, '-o', 'anno_counts', '-g', 'transcript_id', '-T', threads]
    for line in open(sample_list):
          line = line.strip('\n')
          cmd.append(line)
    subprocess.call(cmd)

    # process counts
    template_df = pd.DataFrame(columns = ['Length', 'count', 'factor', 'norm_count'])
    anno_counts = pd.read_csv('anno_counts', sep = '\t', comment = '#', index_col = 0)
    all_genes = pd.unique(all_match['gene'])
    for line in open(sample_list):
        sample = line.strip('\n')
        each_sample_df = pd.DataFrame(columns = ['gene', 'info'])
        each_sample_df_clean = pd.DataFrame(columns = ['gene', 'PSI'])
        print('Calculating expression percentages for: ' + sample)

        # first process anno_counts results
        anno_sub = template_df.copy()
        anno_sub['Length'] = anno_counts['Length']
        anno_sub['count'] = anno_counts[sample]
        # get the median of expressed transcripts
        med = median([i for i in anno_sub['count'].tolist() if i != 0])
        anno_sub['factor'] = 1000 / (anno_sub['Length'] * med)
        # normalize counts for length and median, aks factor
        anno_sub['norm_count'] = anno_sub['count'] * anno_sub['factor']
        anno_sub.set_index(anno_counts.index)

        # then process the RI counts
        RI_counts = sum_table[sum_table['quant_type'] == 'retained_intron']
        RI_sub = template_df.copy()
        RI_sub['count'] = RI_counts[sample]
        RI_sub.set_index(RI_counts.index)
        # the intron length info is in RI, the original count file
        for trans in RI_sub.index:
            RI_sub.at[trans, 'Length'] = RI.at[trans, 'Length']
        # get the median of expressed transcripts
        med = median([i for i in RI_sub['count'].tolist() if i != 0])
        RI_sub['factor'] = 1000 / (RI_sub['Length'] * med)
        # normalize counts for length and median, aka factor
        RI_sub['norm_count'] = RI_sub['count'] * RI_sub['factor']

        # finally the UJ counts
        UJ_counts = sum_table[sum_table['quant_type'] == 'unique_junction']
        UJ_sub = template_df.copy()
        # no need for length correct for UJ
        UJ_sub['count'] = UJ_counts[sample]
        UJ_sub.set_index(UJ_counts.index)
        UJ_sub['Length'] = 1
        # get the median of expressed transcripts
        med = median([i for i in UJ_sub['count'].tolist() if i != 0])
        UJ_sub['factor'] = 1000 / med
        # normalize counts for length and median, aks factor
        UJ_sub['norm_count'] = UJ_sub['count'] * UJ_sub['factor']

        # merge all results for this sample
        sub = pd.concat([anno_sub, RI_sub, UJ_sub], axis = 0)
        # merge to get gene names
        sub = sub.join(all_match, how = 'left')
        # get only expressed transcripts
        sub = sub[sub['count'] != 0]

        # finally process by gene to calculate expression percentage
        for gene in all_genes:
            tmp = sub[sub['gene'] == gene]
            if len(tmp) != 0:
                all_counts = tmp['norm_count'].sum()
                tmp['PSI'] = tmp['norm_count'] * 100 / all_counts
                # the info column would be count;factor;norm_count;PSI
                tmp['info'] = tmp.apply(lambda x: 
                    str(x['count']) + ';' + str(x['factor']) + ';' + 
                    str(x['norm_count']) + ';' + str(x['PSI']), axis = 1)
                # append each gene to each sample df
                each_sample_df = each_sample_df.append(tmp[['gene', 'info']])
                each_sample_df_clean = each_sample_df_clean.append(tmp[['gene', 'PSI']])

        # append each sample to the final tables
        each_sample_df.rename(columns = {'info': sample}, inplace = True)
        each_sample_df_clean.rename(columns = {'PSI': sample}, inplace = True)
        try:
            psi_table = psi_table.join(each_sample_df[sample], how = 'outer')
            psi_table_clean = psi_table_clean.join(each_sample_df_clean[sample], how = 'outer')
        # the first time
        except NameError:
            psi_table = all_match.join(each_sample_df[sample], how = 'left')
            psi_table_clean = all_match.join(each_sample_df_clean[sample], how = 'left')

    psi_table.fillna(0, inplace = True)
    psi_table.sort_values(by = 'gene', inplace = True)  
    psi_table.to_csv('exp_percent_detailed.txt', index = True, index_label = 'transcript', sep = '\t')
    psi_table_clean.fillna(0, inplace = True)
    psi_table_clean.sort_values(by = 'gene', inplace = True)  
    psi_table_clean.to_csv('exp_percent.txt', index = True, index_label = 'transcript', sep = '\t')