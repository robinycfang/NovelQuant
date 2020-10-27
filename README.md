# NovelQuant

## Step 1: find retained introns

`python NovelQuant.py findRI -a annotated.gtf -n novel.gtf`

| Parameter | Description |
|-----------|-------------|
| -a | The gtf file containing information of annotated transcripts. e.g., Gencode annotation |
| -n | The gtf file containing information of novel transcripts. Must have exon information. |

## Step 2: quantify on retained introns

`python NovelQuant.py quantRI -r retained_introns.gtf -l sample_list.txt -p featureCounts_path -t threads`

| Parameter | Description |
|-----------|-------------|
| -r | The gtf file of introns that are retained in novel transcripts. i.e., the output of findRI, retained_introns.gtf |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -p | Path to featureCounts if not in the environment variables |
| -t | Threads to use in featureCounts. Default: 1 |

## Step 3: find unique junctions

`python NovelQuant findUJ -a annotated.gtf -n novel.gtf`

| Parameter | Description |
|-----------|-------------|
| -a | The gtf file containing information of annotated transcripts. e.g., Gencode annotation |
| -n | The gtf file containing information of novel transcripts. Must have exon information. |

## Step 4: quantify on unique junctions

`python NovelQuant quantUJ -a annotated.gtf -n novel.gtf -e uniq_eej.gtf -l sample_list.txt -p featureCounts_path -t threads`

| Parameter | Description |
|-----------|-------------|
| -a | The gtf file containing information of annotated transcripts. e.g., Gencode annotation |
| -n | The gtf file containing information of novel transcripts. Must have exon information. |
| -e | The gtf file of unique junctions. i.e., the output of findUJ, uniq_eej.gtf |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -p | Path to featureCounts if not in the environment variables |
| -t | Threads to use in featureCounts. Default: 1 |

## Step 5: summarize everything

`python NovelQuant sum -r RI_counts.txt -u UJ_counts.txt -l sample_list.txt -p samtools_path -t threads`

| Parameter | Description |
|-----------|-------------|
| -r | Output from quantRI, RI_counts.txt. |
| -u | Output from quantUJ, UJ_counts.txt. |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -p | Path to samtools if not in the environment variables |
| -t | Threads to use in samtools. Default: 1 |

The first column in the final output, NovelQuant_final.txt, represents the names of novel transcripts. The rest of columns represent expression of novel transcripts normalized to total sequencing depth in the samples processed.
