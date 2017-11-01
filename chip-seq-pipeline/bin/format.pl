#!/usr/bin/perl

# grep . AY1278/AY1278.*.qc *.qc | ./format.pl > qc.txt

while (<>) { s/:/ /; @data=split(); shift(@data);
  if (/\s(\S+).filt.nodup.sample.15.SE.tagAlign.gz/) { ($id)=($1); shift(@data); 
    $AY{"$id.tagAlign"}=join(":", @data);
    $AY{"$id.frag_len"}=$data[1];
    $AY{"$id.NSC"}=$data[7];
    $AY{"$id.RSC"}=$data[8];
    $ID{$id}++; } #print "0 $id\n"; }
  elsif (/(\S+).filt.nodup.srt.flagstat.qc (\d+) \+ 0 in total/) { ($id, $total)=($1, $2); $AY{"$id.distinct"}=$total; $ID{$id}++; } #print "1 $id\n"; }
  elsif (/(\S+).filt.nodup.srt.flagstat.qc (\d+) \+ 0 mapped/) { ($id, $mapped)=($1, $2); $AY{"$id.uniqmapped"}=$mapped; $ID{$id}++; } #print "2 $id\n"; }
  elsif (/(\S+).filt.nodup.srt.pbc.qc/) { $pbc=join("|", @data); 
    $AY{"$id.pbc"}=$pbc;
    $AY{"$id.PBC2"}=pop(@data);
    $AY{"$id.PBC1"}=pop(@data);
    $AY{"$id.NRF"}=pop(@data);
    $ID{$id}++; } # print "3 $id\n"; }
  elsif (/(\S+).filt.srt.dup.qc Unknown Library/) { shift(@data); shift(@data); $qc=join(":", @data); $AY{"$id.qc"}=$qc; $ID{$id}++; } #print "4 $id\n"; }
  elsif (/(\S+).raw.srt.flagstat.qc (\d+) \+ (\d+) in total/) { ($id, $rawtotal, $low)=($1,$2,$3); $AY{"$id.hiq_reads"}=$rawtotal; $AY{"$id.loq_reads"}=$low; $ID{$id}++; } #print "5 $id\n"; }
  elsif (/(\S+).raw.srt.flagstat.qc (\d+) \+ 0 mapped/) { ($id, $rawmapped)=($1,$2); $AY{"$id.mapped"}=$rawmapped; $ID{$id}++; } # print "6 $id\n"; }
  elsif (/(\S+).filt.srt.predupmark.qc (\d+) \+ 0 mapped/) {
    ($id, $mapq_filter)=($1,$2); $AY{"$id.mapq_filter"}=$mapq_filter; $ID{$id}++;
  }
}

@tags=("rawtotal", "rawmapped", "uniqtotal", "uniqmapped", "pbc", "qc");
@tags=("hiq_reads", "mapq_filter", "mapped", "fract_mapped", "distinct", "fract_distinct", "NRF", "PBC1", "PBC2", "frag_len", "NSC", "RSC");
@tags=("hiq_reads", "mapped", "mapq_filter", "fract_mapped", "distinct", "fract_distinct", "NRF", "PBC1", "PBC2", "frag_len", "NSC", "RSC");
print "bam"; foreach $tag (@tags) { print "\t$tag"; } print "\n";

#@tags=("NRF", "PBC1", "PBC2", "frag_len", "NSC", "RSC");
#ENCFF245JJL	0.71	0.70	3.07
#pbc=16980381|12112940|8432629|2744414|0.713349|0.696167|3.072652
#qc=16980381:0:0:5292846:0:0:0.311704
#
#	160	1.06	2.48
#AY1278.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc:AY1278.filt.nodup.sample.15.SE.tagAlign.gz	11687535	160	0.168375511858	35	0.1630963	1500	0.1595182	1.055525	2.475443	2

foreach $id (keys %ID) { 
    $key="$id.hiq_reads"; die "::id=$id [$key]?" unless exists $AY{$key};
    $key="$id.mapped"; die "::id=$id [$key]?" unless exists $AY{$key};
    $key="$id.distinct"; die "::id=$id [$key]?" unless exists $AY{$key};
    $key="$id.mapq_filter"; die "::id=$id [$key]?" unless exists $AY{$key};
    next unless exists $AY{"$id.hiq_reads"};
    next unless $AY{"$id.hiq_reads"};
    $AY{"$id.fract_mapped"}=$AY{"$id.mapped"}/$AY{"$id.hiq_reads"}; 
    $AY{"$id.fract_distinct"}=$AY{"$id.distinct"}/$AY{"$id.hiq_reads"}; 
}

foreach $id (keys %ID) {
  print "$id";
  foreach $tag (@tags) { $key="$id.$tag";
    die "::id=$id [$key]?" unless exists $AY{$key};
    print "\t$AY{$key}";
  } print "\n";
}

__END__
AY1278/AY1278	rawtotal=27739210	rawmapped=23516799	uniqtotal=11687535	uniqmapped=11687535	pbc=16980381|12112940|8432629|2744414|0.713349|0.696167|3.072652	qc=16980381:0:0:5292846:0:0:0.311704
bam	hiq_reads	loq_reads	mapped	fract_mapped	distinct	fract_distinct	NRF	PBC1	PBC2	frag_len	NSC	RSC	Notes
ENCFF245JJL	27.7 M	0	23.5 M	0.848	11.7 M	0.421	0.71	0.70	3.07	160	1.06	2.48
AY1278/AY1278	rawtotal=27739210	rawmapped=23516799	uniqtotal=11687535	uniqmapped=11687535	pbc=16980381|12112940|8432629|2744414|0.713349|0.696167|3.072652	qc=16980381:0:0:5292846:0:0:0.311704

AY1278.filt.nodup.srt.flagstat.qc:11687535 + 0 in total (QC-passed reads + QC-failed reads)
AY1278.filt.nodup.srt.flagstat.qc:11687535 + 0 mapped (100.00%:-nan%)
AY1278.filt.nodup.srt.pbc.qc:16980381	12112940	8432629	2744414	0.713349	0.696167	3.072652
AY1278.filt.srt.dup.qc:# picard.sam.markduplicates.MarkDuplicates INPUT=[AY1278.filt.srt.bam] OUTPUT=AY1278.filt.srt.dupmark.bam METRICS_FILE=AY1278.filt.srt.dup.qc REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).* OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
AY1278.filt.srt.dup.qc:LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
AY1278.filt.srt.dup.qc:Unknown Library	16980381	0	0	5292846	0	0	0.311704	
AY1278.raw.srt.flagstat.qc:27739210 + 0 in total (QC-passed reads + QC-failed reads)
AY1278.raw.srt.flagstat.qc:23516799 + 0 mapped (84.78%:-nan%)
