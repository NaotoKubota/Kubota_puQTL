#!/bin/bash

# deepTools
docker run --rm -v /:/wakasugi quay.io/biocontainers/deeptools:3.5.1--py_0 bash -c \
"
    computeMatrix reference-point \
    -S \
    /rhome/naotok/bigdata/ENCODE/GM12878/ENCFF927KAJ_H3K4me3.bigWig \
    /rhome/naotok/bigdata/ENCODE/GM12878/ENCFF564KBE_H3K4me1.bigWig \
    /rhome/naotok/bigdata/ENCODE/GM12878/ENCFF469WVA_H3K27ac.bigWig \
    /rhome/naotok/bigdata/ENCODE/GM12878/ENCFF028KBY_H3K9ac.bigWig \
    /rhome/naotok/bigdata/ENCODE/GM12878/ENCFF919DOR_H3K27me3.bigWig \
    /rhome/naotok/bigdata/ENCODE/GM12878/ENCFF603BJO_ATAC.bigWig \
    -R \
    /rhome/naotok/TSSchoiceQTL/revision/data/puQTL_overlap_eQTL.bed \
    /rhome/naotok/TSSchoiceQTL/revision/data/puQTL_not_eQTL.bed \
    --referencePoint center \
    --skipZeros \
    -b 100000 -a 100000 \
    -p 100 \
    --outFileName /rhome/naotok/bigdata/puQTL/puQTL_overlap_eQTL_matrix.mat.gz && \
    plotHeatmap \
    -m /rhome/naotok/bigdata/puQTL/puQTL_overlap_eQTL_matrix.mat.gz \
    -out /rhome/naotok/TSSchoiceQTL/revision/figure/puQTL_overlap_eQTL_heatmap.png \
    --colorMap viridis \
    --colorList green,blue \
    --heatmapHeight 10 \
    --refPointLabel puQTL \
    --samplesLabel H3K4me3 H3K4me1 H3K27ac H3K9ac H3K27me3 ATAC-seq \
    --regionsLabel puQTL puQTL_eQTL \
    --zMax 4.0 \
    --legendLocation none

"