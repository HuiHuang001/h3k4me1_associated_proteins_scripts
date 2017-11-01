#!/usr/bin/env bash 
#"bin/bwa_aln_fmt_auto_detect.sh"#
set -e
threads=$1
BWA_INDEX_NAME=$2
fastq=$3
sai=$4
log=$5
if [ ${fastq: -3} == ".gz" ]; then 
  input="zcat $fastq"
  elif [ ${fastq: -4} == ".bz2" ]; then 
  input="bzcat $fastq"
  elif [ ${fastq: -6} == ".fastq" ]; then
  input="cat $fastq"
  else 
    echo "File extension of $fastq not recognized. 
          Only *.gz *.bz2 or *.fastq are allowed" 
    exit 1
fi

bwa aln -q 5 -l 32 -k 2 -t $threads $BWA_INDEX_NAME <($input) > $sai 2> $log
bwa samse $BWA_INDEX_NAME $sai <($input) 2>> $log
