#!/bin/bash

##################### Must run in docker container ########################
# docker run --rm -it -v /:/wakasugi naotokubota/qtltools:1.3.1 bash
###########################################################################

for i in 0 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50;
do

	# make covariates
	# sex, top 5 genetic PCs, peer factors
	cat /wakasugi/home/naoto/TSSchoiceQTL/data/txt/E-GEUV-1.sdrf.sex.438.txt \
	<(sed -E 's/ /\t/g' /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/pca/SV.SNP.sort.pca | \
	head -n 6 | sed -e '1d') \
	<(head -n ${i} /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/txt/TMM_peerfactors.txt) | \
	gzip -c \
	> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/covariates_PEER${i}.txt.gz

	# nominal, 1.0 x 10-5
	for j in $(seq 1 22);
	do

		QTLtools cis \
		--vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz \
		--bed /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed.gz \
		--cov /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/covariates_PEER${i}.txt.gz \
		--normal \
		--nominal 0.00001 \
		--chunk ${j} 22 \
		--out /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_PEER${i}_000001_${j}_22.txt &

	done

	wait

	cat /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_PEER${i}_000001_*_22.txt \
	> /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_PEER${i}_000001_full.txt && \
	rm -rf /wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/nominal_PEER${i}_000001_*_22.txt

done