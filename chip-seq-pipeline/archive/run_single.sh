#!/bin/bash

BIN=$(dirname $0)
if [ -z "$ref" ]; then
ref=mm10
fi

# setting up the environment 
java=/opt/jdk1.8.0_72/bin/java
samtools=/usr/bin/samtools
java=/etc/alternatives/jre_1.8.0/bin/java

#
if [ -z "$email" ]; then
email=shz254@ucsd.edu
fi

echo $ref $email

if [ $# -eq 0 ]; then
  exit
fi

function vars {
  if [ -z "$FASTQ_GZ_FILE_1" ]; then echo FASTQ_GZ_FILE_1?; exit; fi
  OFPREFIX=`basename $FASTQ_GZ_FILE_1`
  OFPREFIX=${OFPREFIX%%.fastq.bz2}
  NTHREADS=6
  BWA_INDEX_NAME=/mnt/silencer2/share/bwa_indices/$ref.fa
  SAI_FILE_1="${OFPREFIX}.sai"
  RAW_BAM_PREFIX="${OFPREFIX}.raw.srt"
  RAW_BAM_FILE="${RAW_BAM_PREFIX}.bam" # To be stored
  RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" # QC File
  FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
  FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
  MAPQ_THRESH=30
  TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.predupmark.bam"
  MARKDUP="/srv/gs1/software/picard-tools/1.92/MarkDuplicates.jar"
  MARKDUP="/opt/picard-tools-1.129/picard.jar MarkDuplicates"
  MARKDUP="/mnt/silencer2/share/picard-tools-2.1.0/picard.jar MarkDuplicates"
  DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file
  FINAL_BAM_PREFIX="${OFPREFIX}.filt.nodup.srt"
  FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
  FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
  FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file
  PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"
  FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"
  NREADS=15000000
  SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.$((NREADS / 1000000)).SE.tagAlign.gz"
  CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
  CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"

  if [ -z "$OFPREFIX" ]; then echo OFPREFIX?; exit; fi
  # compgen -A variable 
  # compgen -A variable | awk '{if(n)print}/^_$/{n++}'
}

# map fastq to the genome, and filter reads. 
for FASTQ_GZ_FILE_1 in $*; do vars;
  if [ -s ${RAW_BAM_FILE_MAPSTATS} ]; then
    echo "pass 1 ${RAW_BAM_FILE_MAPSTATS}"
    continue
  fi
  echo ${RAW_BAM_FILE_MAPSTATS}
  echo "$OFPREFIX::"; set | grep $OFPREFIX; echo
  touch ${RAW_BAM_FILE_MAPSTATS}

  bwa aln -q 5 -l 32 -k 2 -t ${NTHREADS} ${BWA_INDEX_NAME} <(bzcat ${FASTQ_GZ_FILE_1}) > ${SAI_FILE_1}
  bwa samse ${BWA_INDEX_NAME} ${SAI_FILE_1} <(bzcat ${FASTQ_GZ_FILE_1}) | $samtools view -Su - | $samtools sort - ${RAW_BAM_PREFIX}
  rm ${SAI_FILE_1}
  $samtools flagstat ${RAW_BAM_FILE} > ${RAW_BAM_FILE_MAPSTATS}
  
  $samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} > ${TMP_FILT_BAM_FILE}
  $samtools flagstat ${TMP_FILT_BAM_FILE} > ${TMP_FILT_BAM_FILE%%bam}qc
done

### mark duplicates. 
for FASTQ_GZ_FILE_1 in $*; do vars
  if [ -s ${DUP_FILE_QC} ]; then
    echo "pass 2 ${DUP_FILE_QC}"
    continue
  fi
  echo ${DUP_FILE_QC}
  echo "$OFPREFIX::"; set | grep $OFPREFIX; echo
  touch ${DUP_FILE_QC}
 
   $java -Xmx4G -jar ${MARKDUP} TMP_DIR=`pwd`/tmp INPUT=${TMP_FILT_BAM_FILE} OUTPUT=${FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false 
done

## remove duplicates and produce final bam file. 
for FASTQ_GZ_FILE_1 in $*; do vars
  if [ -s ${PBC_FILE_QC} ]; then
    echo "pass 3 ${PBC_FILE_QC}"
    continue
  fi
  echo ${PBC_FILE_QC}
  echo "$OFPREFIX::"; set | grep $OFPREFIX; echo
  touch ${PBC_FILE_QC}

  $samtools view -F 1804 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}
  # Index Final BAM file
  $samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}
  $samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}
  bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort -T . | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
  # rm ${FILT_BAM_FILE}
done


export PATH=/mnt/silencer2/share/phantompeakqualtools/:/usr/local/bin/:/mnt/silencer2/share/ENCODE/bin:$PATH

## calculate phantom score. 

for FASTQ_GZ_FILE_1 in $*; do vars
  if [ -e ${CC_SCORES_FILE}.flag ]; then
    echo "pass 4 ${CC_SCORES_FILE}"
    continue
  fi
  echo ${CC_SCORES_FILE}
  echo "$OFPREFIX::"; set | grep $OFPREFIX; echo
  echo $FASTQ_GZ_FILE_1 > ${CC_SCORES_FILE}.flag

  if ! [ -e ${FINAL_TA_FILE} ]; then
    bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${FINAL_TA_FILE}
  fi
  if ! [ -e  ${SUBSAMPLED_TA_FILE} ]; then
    zcat ${FINAL_TA_FILE} | grep -v "chrM" | shuf -n ${NREADS} | gzip -c > ${SUBSAMPLED_TA_FILE}
  fi
  # CC_SCORE FILE format
  # Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
  echo "Rscript /mnt/silencer2/share/phantompeakqualtools/run_spp_nodups.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE} -tmpdir=."
  Rscript /mnt/silencer2/share/phantompeakqualtools/run_spp_nodups.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE} -tmpdir=.
  sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
  mv temp ${CC_SCORES_FILE}
done

#rm ${TMP_FILT_BAM_FILE} ${RAW_BAM_FILE} ${FILT_BAM_FILE} ${FINAL_BAM_FILE}

## convert bam files to bigwig files. 
m=`ls *.filt.nodup.srt.bam 2>/dev/null|wc -l`
if [ $m -ne 0 ]; then
for FASTQ_GZ_FILE_1 in $*; do vars
for bam in $OFPREFIX*.filt.nodup.srt.bam; do
  id=`basename $bam`
  id=${id%%.bam}
  ext_len=`$samtools view $bam | awk '{print $10; exit}' | awk '{print 300-length}'`
  gs=$id.gs
  bw=$id.bw
  if ! [ -e $bw ]; then
    echo $id $ext_len
    $samtools view -H $bam | awk '$1 == "@SQ" {OFS="\t";print $2,$3}' - | sed 's/.N://g' > $gs;
    $samtools view -b $bam | bedtools bamtobed | bedtools slop -s -l 0 -r $ext_len -i stdin -g $gs |
      bedtools genomecov -g $gs -i stdin -bg | wigToBigWig stdin $gs $bw; rm $gs
  else
    echo Find $bw, skip
  fi
done

done
fi

#rm *.fastq.gz.filt.srt.*bam
# rm ${TMP_FILT_BAM_FILE} ${RAW_BAM_FILE} ${FILT_BAM_FILE} ${FINAL_BAM_FILE}

## last step

(echo | ../format.pl
for FASTQ_GZ_FILE_1 in $*; do vars
  grep . $OFPREFIX*.qc | $BIN/format.pl | grep -v fract_mapped
done) > qc.txt

echo "done ${OFPREFIX}"

#if [[ -t 1 ]]; then
#  $BIN/format.sh
#else
#  $BIN/format.sh | mail -s "QC done" $email
#fi
