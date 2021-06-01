# NovelQuant

NovelQuant takes GTF files of annotated and novel transcripts as inputs, and identifies unique regions (retained introns and unique junctions) in the novel transcripts, relative to annotated transcripts. NovelQuant then works on BAM files (can be either short and long reads) of the samples of your interest and counts reads falling upon the unique regions to quantify the novel transcripts.

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

`python NovelQuant sum -r RI_counts.txt -u UJ_counts.txt -l sample_list.txt -p samtools_path`

| Parameter | Description |
|-----------|-------------|
| -r | Output from quantRI, RI_counts.txt. |
| -u | Output from quantUJ, UJ_counts.txt. |
| -l | A list of BAM file(s) to be processed. Each line should be the path of each BAM file. |
| -p | Path to samtools if not in the environment variables |

The first column in the final output, NovelQuant_final.txt, represents the names of novel transcripts. The rest of columns represent expression of novel transcripts normalized to total sequencing depth in the samples processed.
