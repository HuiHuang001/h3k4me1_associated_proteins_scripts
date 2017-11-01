#!/bin/bash

#grep . AY1278/AY1278.*.qc | ./format.pl > qc.txt
dir=`dirname $0`
#if ! [ -e qc.txt ]; then
#grep . *.qc | $dir/format.pl > qc.txt
#fi

path=`pwd | sed 's/.home.bli.//'`
me=`whoami`
md=`date +"%m%d"`
dir=/mnt/silencer2/home/$me/public_html/ENCODE/`date +"%m%d"`
if ! [ -e "$dir" ]; then
  mkdir -p $dir
fi
rsync -av *.bw $dir/.

for bw in *.bw; do
  echo "track type=bigWig name=\"$bw\" description=\"$bw\" bigDataUrl=http://enhancer.sdsc.edu/$me/ENCODE/$md/$bw visibility=full"
done > tracks.txt

echo "<pre>"
cat qc.txt 
echo
cat tracks.txt
echo "</pre>"
