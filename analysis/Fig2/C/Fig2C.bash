#!/bin/bash

################### Must run in docker container #####################
# docker run --rm -it -v /:/wakasugi naotokubota/qtltools:1.3.1 bash
######################################################################

### eQTL ###

# 1. make input files

# significant eQTL (best hit)
zcat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_full.txt.gz | \
awk -F" " -v OFS="\t" '$23 == 1 && $12 == 0{print $9,$10-1,$11,$8,$1,$5}' | \
sort -k1,1 -k2,2n \
> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/bed/conditional_PEER50_01_besthit.bed
# all phenotypes
zcat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed.gz | \
sed -e 1d | \
awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$4,".",$6}' \
> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/bed/TMM_simple.bed

# 2. fenrich

# Histone and ATAC
rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_histone.txt

for MARK in H3K4me3 H3K4me1 H3K27ac H3K9ac H3K27me3 ATAC
do

	QTLtools fenrich \
	--qtl /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/bed/conditional_PEER50_01_besthit.bed \
	--tss /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/bed/TMM_simple.bed \
	--bed /wakasugi/mnt/data5/naoto/ENCODE/GM12878/${MARK}.sort.bed.gz \
	--permute 1000 \
	--out /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_${MARK}.txt && \
	cat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_${MARK}.txt | \
	awk -F" " -v OFS=" " -v MARK=${MARK} '{print $0,MARK}' | \
	sed -e 's/ /\t/g' \
	>> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_histone.txt && \
	rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_${MARK}.txt

done

# TF footprints
rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_TFfootprints.txt

for ID in GM06990-DS7748 GM12865-DS12436
do

	CELL=`echo ${ID} | cut -d"-" -f 1`

	QTLtools fenrich \
	--qtl /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/bed/conditional_PEER50_01_besthit.bed \
	--tss /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/bed/TMM_simple.bed \
	--bed /wakasugi/mnt/data5/naoto/TFfootprints/${ID}/interval.all.fps.0.05.bed \
	--permute 1000 \
	--out /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_${CELL}.txt && \
	cat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_${CELL}.txt | \
	awk -F" " -v OFS=" " -v CELL=${CELL} '{print $0,CELL}' | \
	sed -e 's/ /\t/g' \
	>> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_TFfootprints.txt && \
	rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_${CELL}.txt

done

cat \
/wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_histone.txt \
/wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_TFfootprints.txt \
> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment.txt
