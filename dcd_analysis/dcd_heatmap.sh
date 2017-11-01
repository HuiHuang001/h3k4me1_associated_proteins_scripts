#!/bin/bash

computeMatrix reference-point --referencePoint center -S ../sig_FE_pool/r1_me1.pool.FE.bw ../sig_FE_pool/dcd_me1.pool.FE.bw -R me1.dist.decreased.bed me1.dist.unchanged.bed -a 10000 -b 10000 -p 24 -o me1.dist.new.matrix.gz

plotHeatmap -m me1.dist.new.matrix.gz --sortUsingSamples 1 --colorList 'white,white,red,red' --whatToShow 'heatmap and colorbar' -out me1.dist.new.heatmap.pdf --outFileSortedRegions me1.dist.heatmap.new.sorted.bed --legendLocation "lower-center" --refPointLabel "center"

plotProfile -m me1.dist.new.matrix.gz --perGroup -out me1.dist.new.profile.pdf --colors black red -–numPlotsPerRow 1 --refPointLabel "center"

for name in k27ac brg1 dpf2 me2 me3
do

computeMatrix reference-point --referencePoint center -S ../sig_FE_pool/r1_$name.pool.FE.bw ../sig_FE_pool/dcd_$name.pool.FE.bw -R me1.dist.heatmap.new.sorted.bed --sortRegions keep -a 10000 -b 10000 -p 24 -o $name.dist.new.matrix.gz

plotHeatmap -m $name.dist.new.matrix.gz --sortRegions no --colorList 'white,white,red,red' --whatToShow 'heatmap and colorbar' -out $name.dist.new.heatmap.pdf --legendLocation "lower-center" --refPointLabel "center"

plotProfile -m $name.dist.new.matrix.gz --perGroup -out $name.dist.new.profile.pdf --colors black red --numPlotsPerRow 1 --refPointLabel "center"

done

