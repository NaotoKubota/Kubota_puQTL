#!/bin/bash

# make puQTL bed
awk -F" " -v OFS="\t" '$23 == 1{print $9,$10-1,$10,$8}' \
<(zcat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz) | \
sort -k1,1 -k2,2n \
> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/bed/conditional_PEER25_01_besthit.bed

# deepTools
docker run --rm -v /:/wakasugi quay.io/biocontainers/deeptools:3.5.1--py_0 bash -c "\
computeMatrix reference-point \
-S \
/wakasugi/mnt/data5/naoto/ENCODE/GM12878/ENCFF927KAJ_H3K4me3.bigWig \
/wakasugi/mnt/data5/naoto/ENCODE/GM12878/ENCFF564KBE_H3K4me1.bigWig \
/wakasugi/mnt/data5/naoto/ENCODE/GM12878/ENCFF469WVA_H3K27ac.bigWig \
/wakasugi/mnt/data5/naoto/ENCODE/GM12878/ENCFF028KBY_H3K9ac.bigWig \
/wakasugi/mnt/data5/naoto/ENCODE/GM12878/ENCFF919DOR_H3K27me3.bigWig \
/wakasugi/mnt/data5/naoto/ENCODE/GM12878/ENCFF603BJO_ATAC.bigWig \
-R /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/bed/conditional_PEER25_01_besthit.bed \
--referencePoint center \
--skipZeros \
-b 100000 -a 100000 \
-p 24 \
--outFileName /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/tsv/conditional_PEER25_01_SNP_GM12878_histone_matrix.mat.gz\
"

docker run --rm -v /:/wakasugi quay.io/biocontainers/deeptools:3.5.1--py_0 bash -c "\
plotHeatmap \
-m /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/tsv/conditional_PEER25_01_SNP_GM12878_histone_matrix.mat.gz \
-out /wakasugi/home/naoto/TSSchoiceQTL/png_svg/conditional_PEER25_01_SNP_GM12878_histone_heatmap.png \
--colorMap viridis \
--heatmapHeight 10 \
--refPointLabel puQTL \
--samplesLabel H3K4me3 H3K4me1 H3K27ac H3K9ac H3K27me3 ATAC-seq \
--zMax 1.0 \
--legendLocation none\
"