#!/bin/bash

for name in brg1 chd1 dpf2 p300 me1 me2 me3 smc3 k27ac

do

###define reproducible peaks as peaks that were found in pool replicates as well as in each individual replicate.



bedtools intersect -a ./peaks_pooled/dcd_$name.pool_peaks.narrowPeak -b ./peaks/dcd_$name.rep1_peaks.narrowPeak -u > tmp1.bed

bedtools intersect -a ./peaks_pooled/dcd_$name.pool_peaks.narrowPeak -b ./peaks/dcd_$name.rep2_peaks.narrowPeak -u > tmp2.bed

bedtools intersect -a tmp1.bed -b tmp2.bed > ./repro_peaks/dcd_$name.repro.peak

rm tmp1.bed
rm tmp2.bed



bedtools intersect -a ./peaks_pooled/r1_$name.pool_peaks.narrowPeak -b ./peaks/r1_$name.rep1_peaks.narrowPeak -u > tmp1.bed

bedtools intersect -a ./peaks_pooled/r1_$name.pool_peaks.narrowPeak -b ./peaks/r1_$name.rep2_peaks.narrowPeak -u > tmp2.bed

bedtools intersect -a tmp1.bed -b tmp2.bed > ./repro_peaks/r1_$name.repro.peak

rm tmp1.bed
rm tmp2.bed
### merge peaks in dcd and wt (r1) for comparative analysis

cat ./repro_peaks/dcd_$name.repro.peak ./repro_peaks/r1_$name.repro.peak | sort -k1,1 -k2,2n > ./repro_peaks/merge.tmp.bed

bedtools merge -i ./repro_peaks/merge.tmp.bed -d 50 | sort -k1,1 -k2,2n | awk -v var="$name" 'BEGIN{OFS="\t";i=1}{print $1,$2,$3,var"_merged_peak"i;i++}' > ./repro_peaks/$name.merge.peak 

rm merge.tmp.bed

### make saf file
cat ./repro_peaks/$name.merge.peak | awk 'BEGIN{OFS="\t";print "GeneID\tChr\tStart\tEnd\tStrand"}{print $1"."$2"."$3,$1,$2,$3,"."}' > ./repro_peaks/$name.merge.saf

### count reads fall in the repro peaks for each r1, dcd replicate

featureCounts -a ./repro_peaks/$name.merge.saf -F SAF -T 10 -o ./repro_peaks/$name.count.txt ./bam/r1_$name.rep1.filt.nodup.srt.bam ./bam/r1_$name.rep2.filt.nodup.srt.bam ./bam/dcd_$name.rep1.filt.nodup.srt.bam ./bam/dcd_$name.rep2.filt.nodup.srt.bam 

### generate signal score

bigWigAverageOverBed ./sig_FE_pool/dcd_$name.pool.FE.bw ./repro_peaks/$name.merge.peak ./repro_peaks/dcd_$name.score.txt

#bigWigAverageOverBed ./sig_FE_pool/dko_$name.pool.FE.bw ./repro_peaks/$name.merge.peak ./repro_peaks/dko_$name.score.txt

bigWigAverageOverBed ./sig_FE_pool/r1_$name.pool.FE.bw ./repro_peaks/$name.merge.peak ./repro_peaks/r1_$name.score.txt

done

###find all gencode annotated tss

zcat gencode.vM1.annotation.gtf.gz | awk 'BEGIN{OFS="\t"}{if($3=="transcript" && $7=="+")print $1,$4,$4+1;if($3=="transcript" && $7=="-")print $1,$5,$5+1}' | sort -k1,1 -k2,2n > mm9.all.tss.bed

### if the closet tss is 2kb away from a peak, that peak is defined as distal peak and vice versa.

bedtools closest -a ./repro_peaks/me1.merge.peak -b mm9.all.tss.bed | awk 'BEGIN{OFS="\t"}{if(($3+$2)/2-$6>2000 || ($3+$2)/2-$6<-2000) print $1,$2,$3,$4}' | sort -k1,1 -k2,2n > ./k4me1_region/me1.tmp.dist.peak

bedtools intersect -a ./repro_peaks/me1.merge.peak -b ./k4me1_region/me1.tmp.dist.peak -u | sort -k1,1 -k2,2n > ./k4me1_region/me1.merge.dist.peak

bedtools intersect -a ./repro_peaks/me1.merge.peak -b ./k4me1_region/me1.merge.dist.peak -v | sort -k1,1 -k2,2n > k4me1_region/me1.merge.prox.peak

bedtools intersect -a ./k4me1_region/me1.merge.dist.peak -b ./repro_peaks/k27ac.merge.peak -v | sort -k1,1 -k2,2n > ./k4me1_region/me1.primed.peak

bedtools intersect -a ./k4me1_region/me1.merge.dist.peak -b ./repro_peaks/k27ac.merge.peak -u | sort -k1,1 -k2,2n > ./k4me1_region/me1.active.peak
