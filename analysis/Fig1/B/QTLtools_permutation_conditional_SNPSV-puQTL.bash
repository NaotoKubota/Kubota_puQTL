#!/bin/bash

##################### Must run in docker container ########################
# docker run --rm -it -v /:/wakasugi naotokubota/qtltools:1.3.1 bash
###########################################################################

 # Set PEER number
# PEER=PEER25
PEER=PEER0

# permutation
for j in $(seq 1 22);
do

    QTLtools cis \
    --vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz \
    --bed /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed.gz \
    --cov /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/covariates_${PEER}.txt.gz \
    --normal \
    --std-err \
    --permute 1000 \
    --chunk ${j} 22 \
    --out /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/permutation_${PEER}_${j}_22.txt &

done

wait

cat /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/permutation_${PEER}_*_22.txt | gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/permutation_${PEER}_full.txt.gz && \
rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/permutation_${PEER}_*_22.txt && \
Rscript /qtltools/scripts/qtltools_runFDR_cis.R \
/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/permutation_${PEER}_full.txt.gz 0.1 \
/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/permutation_${PEER}_01

# conditional
for j in $(seq 1 22);
do

    QTLtools cis \
    --vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz \
    --bed /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed.gz \
    --cov /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/covariates_${PEER}.txt.gz \
    --mapping /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/permutation_${PEER}_01.thresholds.txt \
    --normal \
    --std-err \
    --chunk ${j} 22 \
    --out /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_${PEER}_01_${j}_22.txt &

done

wait

cat /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_${PEER}_01_*_22.txt | gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_${PEER}_01_full.txt.gz && \
rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_${PEER}_01_*_22.txt
