#!/bin/bash

###################### Must run in docker container #######################
# docker run --rm -it -v /:/wakasugi naotokubota/qtltools:1.3.1 bash
###########################################################################

# puQTL
zcat /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz | \
awk '$23 == 1' | \
gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_besthit.txt.gz

for j in $(seq 0 22)
do

	QTLtools rtc \
	--vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz \
	--bed /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed.gz \
    --cov /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/covariates_PEER25.txt.gz \
	--hotspots /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh38_cM10_merge.bed \
	--gwas-cis \
	/wakasugi/mnt/data5/naoto/GWASCatalog/gwas_catalog_v1.0.2-associations_e104_r2021-12-07_allSNPs.txt \
	/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_besthit.txt.gz \
	--normal \
	--chunk ${j} 22 \
	--out /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/rtc_all_PEER25_${j}_22.txt &

done

wait

cat /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/rtc_all_PEER25_*_22.txt | gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/rtc_all_PEER25_full.txt.gz && \
rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/rtc_all_PEER25_*_22.txt

# eQTL
zcat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_full.txt.gz | \
awk '$23 == 1' | \
gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_besthit.txt.gz
for j in $(seq 0 22)
do

	QTLtools rtc \
	--vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz \
	--bed /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed.gz \
    --cov /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/covariates_PEER50.txt.gz \
	--hotspots /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh38_cM10_merge.bed \
	--gwas-cis \
	/wakasugi/mnt/data5/naoto/GWASCatalog/gwas_catalog_v1.0.2-associations_e104_r2021-12-07_allSNPs.txt \
	/wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_besthit.txt.gz \
	--normal \
	--chunk ${j} 22 \
	--out /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/rtc_all_PEER50_${j}_22.txt &

done

wait

cat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/rtc_all_PEER50_*_22.txt | gzip -c \
> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/rtc_all_PEER50_full.txt.gz && \
rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/rtc_all_PEER50_*_22.txt
