# NovelQuant

NovelQuant takes GTF files of annotated and novel transcripts as inputs, and identifies unique regions (retained introns and unique junctions) in the novel transcripts, relative to annotated transcripts. NovelQuant then works on BAM files (can be either short and long reads) of the samples of your interest and counts reads falling upon the unique regions to quantify the novel transcripts. Optionally, under a statistical assumption, an expression normalization can be conducted to calculate expression percentages of annotated and novel transcripts within each gene, which enables their relative comparision. 

## Step 1: find retained introns

NovelQuant compares GTF files of annotated and novel transcripts, and searches for intronic regions that are retained in the novel transcripts.

`python NovelQuant.py findRI -a annotated.gtf -n novel.gtf`

| Parameter | Description |
|-----------|-------------|
| -a | The gtf file containing information of annotated transcripts. e.g., Gencode annotation |
| -n | The gtf file containing information of novel transcripts. Must have exon information. |

## Step 2: quantify on retained introns

Uses featureCounts to count reads falling upon retained introns identified in the Step 1.

`python NovelQuant.py quantRI -r retained_introns.gtf -l sample_list.txt -p featureCounts_path -t threads`

| Parameter | Description |
|-----------|-------------|
| -r | The gtf file of introns that are retained in novel transcripts. i.e., the output of findRI, retained_introns.gtf |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -p | Path to featureCounts if not in the environment variables |
| -t | Threads to use in featureCounts. Default: 1 |

## Step 3: find unique junctions

NovelQuant compares GTF files of annotated and novel transcripts, and searches for exon-exon junctions that are unique in the novel transcripts.

`python NovelQuant findUJ -a annotated.gtf -n novel.gtf`

| Parameter | Description |
|-----------|-------------|
| -a | The gtf file containing information of annotated transcripts. e.g., Gencode annotation |
| -n | The gtf file containing information of novel transcripts. Must have exon information. |

## Step 4: quantify on unique junctions

Uses featureCounts to count junction reads spanning the unique junctions that are identified in Step 3.

`python NovelQuant quantUJ -n novel.gtf -e uniq_eej.gtf -l sample_list.txt -p featureCounts_path -t threads`

| Parameter | Description |
|-----------|-------------|
| -n | The gtf file containing information of novel transcripts. Must have exon information. |
| -e | The gtf file of unique junctions. i.e., the output of findUJ, uniq_eej.gtf |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -p | Path to featureCounts if not in the environment variables |
| -t | Threads to use in featureCounts. Default: 1 |

## Step 5: summarize everything

Merges results of quantRI and quantUJ, and normalize to total sequencing depth of each sample.

`python NovelQuant sum -r RI_counts.txt -u UJ_counts.txt -l sample_list.txt -st samtools_path`

| Parameter | Description |
|-----------|-------------|
| -r | Output from quantRI, RI_counts.txt. |
| -u | Output from quantUJ, UJ_counts.txt. |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -st | Path to samtools if not in the environment variables |

The first column in the final output, `NovelQuant_final.txt`, represents the names of novel transcripts. The rest of columns represent expression levels of novel transcripts in each sample (corrected for sequencing depth).

Optionally, follow this command to calculate expression percentages of annotated and novel transcripts in each gene.

`python NovelQuant sum -r RI_counts.txt -u UJ_counts.txt -l sample_list.txt -st samtools_path --CalExpPerc -a annotated.gtf -n novel.gtf -fc featureCounts_path -t threads`

| Parameter | Description |
|-----------|-------------|
| -r | Output from quantRI, RI_counts.txt. |
| -u | Output from quantUJ, UJ_counts.txt. |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -st | Path to samtools if not in the environment variables |
| --CalExpPerc | Calculate expression percentages of annotaed and novel transcripts in each gene |
| -a | The gtf file containing information of annotated transcripts. e.g., Gencode annotation |
| -n | The gtf file containing information of novel transcripts. Must have exon information. |
| -fc | Path to featureCounts if not in the environment variables |
| -t | Threads to use in featureCounts. Default: 1 |

A direct comparison of expression levels between annotated and novel transcripts is not feasible due to a unique quantification method used for novel transcripts by NovelQuant. We perform an assumption that for the expression data, selection of annotated and novel transcripts are two random samplings from the transcriptome, so they follow the same statistical distribution, and the true expression levels of median annotated and novel transcripts are identical. Therefore, NovelQuant uses these median transcripts as a reference to normalize annotated and novel transcripts separately (calcluate relative expression to the reference), which eventually makes the comparison between annotated and novel transcripts possible. This optional step generates 2 outputs. For each column sample in `exp_percent_detailed.txt`, each entry shows original counts, factor, normalized counts and expresison percentage, delimited by comma. For each sample column in `exp_percent.txt`, each entry only shows the final expresison percentage for ease of the downstream analysis.
