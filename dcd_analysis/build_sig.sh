#!/bin/bash

#PBS -q hotel
#PBS -N hui_job
#PBS -l nodes=1:ppn=1,pmem=8gb,walltime=24:00:00
#PBS -V
#PBS -M huh025@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -t 2-27

## load module
module load python
## running programs


chip=$(sed -n ${PBS_ARRAYID}p /oasis/tscc/scratch/huh025/HH_CHIP_MLL/pool.info.txt | cut -f 1)
input=$(sed -n ${PBS_ARRAYID}p /oasis/tscc/scratch/huh025/HH_CHIP_MLL/pool.info.txt | cut -f 2)
name=$(sed -n ${PBS_ARRAYID}p /oasis/tscc/scratch/huh025/HH_CHIP_MLL/pool.info.txt | cut -f 3)
inpath="/oasis/tscc/scratch/huh025/HH_CHIP_MLL/peaks_pooled"
outpath="/oasis/tscc/scratch/huh025/HH_CHIP_MLL/sig_FE_pooled"

### In pool.info.txt is information for pooled replicates.
### eg. r1_me2.pool_treat_pileup.bdg	r1_me2.pool_control_lambda.bdg	r1_me2.pool
if [ ! -d $outpath ]; then mkdir $outpath;
fi
#mkdir $outpath

## running MACS2

echo $name

##build signal file for pooled replicates
macs2 bdgcmp -t $inpath/$chip -c $inpath/$input -m FE -o $name --outdir $outpath
sort $outpath/$name -k1,1 -k2,2n > $outpath/${name}.sorted

bedGraphToBigWig $outpath/${name}.sorted $outpath/mm9.chrom.sizes $outpath/${name}.FE.bw

rm $outpath/$name



