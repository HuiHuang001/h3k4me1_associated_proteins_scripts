#!/bin/bash

#PBS -q hotel
#PBS -N hui_job
#PBS -l nodes=1:ppn=1,pmem=8gb,walltime=24:00:00
#PBS -V
#PBS -M huh025@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -t 38-109

## load module
module load python
## running programs


chip=$(sed -n ${PBS_ARRAYID}p /oasis/tscc/scratch/huh025/HH_CHIP_MLL/sample.info.txt | cut -f 1)
input=$(sed -n ${PBS_ARRAYID}p /oasis/tscc/scratch/huh025/HH_CHIP_MLL/sample.info.txt | cut -f 2)
name=$(sed -n ${PBS_ARRAYID}p /oasis/tscc/scratch/huh025/HH_CHIP_MLL/sample.info.txt | cut -f 3)
inpath="/oasis/tscc/scratch/huh025/HH_CHIP_MLL/bam"
outpath="/oasis/tscc/scratch/huh025/HH_CHIP_MLL/peaks"
###In sample.info.txt are names of sample and mapped chip and input bam files.
### eg. r1_me1.rep2.filt.nodup.srt.bam	r1_input.rep2.filt.nodup.srt.bam	r1_me1.rep2

if [ ! -d $outpath]; then
mkdir $outpath
fi
#mkdir $outpath/$name

## running MACS2 to call peaks with default parameters

macs2 callpeak -t $inpath/$chip -c $inpath/$input -f BAM -g mm -m 5 50 -B -n $name --outdir $outpath

###generate signal track using bamCoverage
bamCoverage -b $inpath/$chip -o $name --normalizeUsingRPKM -e 300 -bs 25 --smoothLength 75
