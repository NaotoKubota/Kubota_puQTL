#!/bin/bash

################### Must run in docker container #####################
# docker run --rm -it -v /:/wakasugi naotokubota/qtltools:1.3.1 bash
######################################################################

rm -rf \
/wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_chromHMM.txt \

for i in {1..25}
do


	### eQTL ###

	QTLtools fenrich \
	--qtl /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/bed/conditional_PEER50_01_besthit.bed \
	--tss /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM_simple.bed \
	--bed /wakasugi/mnt/data5/naoto/chromHMM/E116_GM12878_Lymphoblastoid_Cell_Line_25_imputed12marks_hg38lift_mnemonics_${i}.bed \
	--permute 1000 \
	--out /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_chromHMM_${i}.txt && \
	cat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_chromHMM_${i}.txt | \
	awk -F" " -v OFS=" " -v GROUP=${i} '{print $0,GROUP}' | \
	sed -e 's/ /\t/g' \
	>> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_chromHMM.txt && \
	rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_chromHMM_${i}.txt


done
