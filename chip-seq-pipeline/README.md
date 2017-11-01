# chip-seq-pipeline
This is the ChIP-seq mapping pipeline for our lab. The script is still at the very preliminary stage, please let us know how to improve!

## How to run it.
### Single sample. 
1. Install snakemake 
2. Go to your work directory, create and copy/link your fastq files in the fastq
folder
3. <path-to-the-pipeline>/run_chipseq_vanilla.sh -g hg19 -e [email@gmail.com]. It will process every sample fastqs in the fastq/ directory, and send an e-mail when it is done. 


## Output files:
All the output files will be relative to your project directory. 
1. bam/[name].filt.nodup.srt.bam: the final bam files after filter, duplicate removal, and sort. 
2. bigWig/[name].filt.nodup.srt.bw: bigWig file for the bam.
3. qc: summary statistics for the bam files.
    1. \*.flagstat.qc: summary stats for different steps of bam files.
    2. [name].\*.cc.qc: phantom peak summary.
    3. [name].\*.pbc.qc: PBC summary. 
    4. [name].\*.plot.pdf: tag shift size estimate plot.
4. logs/: running logs for different parts of the pipeline


## Dependencies
* BWA
* samtools
* bedtools
* picard
* phantompeaktools


## Plans for implementation
At its complete form, this pipeline will include work flow for analysis of (1) a single sample, and (2) replicates using the ENCODE IDR method.

## Roadmap
2017.05.26: The Snakemake pipeline for single sample is available.
2017.05.17: Currently we have only the workflow for a single sample. ***And it is not modularized yet***. 
