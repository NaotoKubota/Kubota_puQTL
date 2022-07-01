#!/bin/bash

##################### Must run in docker container ########################
# docker run --rm -it -v /:/wakasugi naotokubota/qtltools:1.3.1 bash
###########################################################################

 # Set PEER number
PEER=PEER50

for j in $(seq 1 22);
do

    QTLtools cis \
    --vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz \
    --bed /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed.gz \
    --cov /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/covariates_${PEER}.txt.gz \
    --normal \
    --nominal 1.0 \
    --std-err \
    --chunk ${j} 22 \
    --out /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_${PEER}_${j}_22.txt &

done

wait

cat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_${PEER}_*_22.txt | gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_${PEER}_full.txt.gz && \
rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_${PEER}_*_22.txt && \
zcat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_${PEER}_full.txt.gz | \
awk -F" " '$2 == "chr22"{print}' | gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_${PEER}_chr22.txt.gz
