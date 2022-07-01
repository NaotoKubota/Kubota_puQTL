#!/bin/bash

### Must run in docker container! ###
# docker run -v /:/wakasugi -it --rm aylab/fithichip bash

BAMPATH=/wakasugi/mnt/data4/naoto/puQTL/HiChIP/HiC-Pro_results_hg38/bowtie_results/bwt2 && \
for SAMPLE in SRR5831489 SRR5831490
do

	mkdir -p /wakasugi/mnt/data4/naoto/puQTL/HiChIP/FitHiChIP/macs2/${SAMPLE} && \
	macs2 callpeak \
	-t ${BAMPATH}/${SAMPLE}/${SAMPLE}_filtered_hg38.bwt2pairs.bam \
	-f BAM \
	-g hs \
	-n ${SAMPLE} \
	-B \
	-q 0.01 \
	--outdir /wakasugi/mnt/data4/naoto/puQTL/HiChIP/FitHiChIP/macs2/${SAMPLE} \
	> /wakasugi/mnt/data4/naoto/puQTL/HiChIP/FitHiChIP/macs2/${SAMPLE}/macs2.log 2>&1 && \
	FITHICHIPDIR=/wakasugi/mnt/data4/naoto/puQTL/HiChIP/FitHiChIP/results/${SAMPLE} && \
	mkdir -p ${FITHICHIPDIR} && \
	sed \
	-e "s/ValidPairs=\.\/TestData\/Sample_ValidPairs\.txt\.gz/ValidPairs=\/wakasugi\/mnt\/data4\/naoto\/puQTL\/HiChIP\/HiC-Pro_results_hg38\/hic_results\/data\/${SAMPLE}\/${SAMPLE}\.allValidPairs/" \
	-e "s/PeakFile=\.\/TestData\/Sample\.Peaks\.gz/PeakFile=\/wakasugi\/mnt\/data4\/naoto\/puQTL\/HiChIP\/FitHiChIP\/macs2\/${SAMPLE}\/${SAMPLE}_peaks\.narrowPeak/" \
	-e "s/OutDir=\.\/TestData\/results\//OutDir=\/wakasugi\/mnt\/data4\/naoto\/puQTL\/HiChIP\/FitHiChIP\/results\/${SAMPLE}\//" \
	-e "s/ChrSizeFile=\.\/TestData\/chrom_hg19\.sizes/ChrSizeFile=\/wakasugi\/mnt\/data6\/naoto\/twins1\/RefGenome\/hg38_genomeSize\.txt/" \
	-e "s/PREFIX=FitHiChIP/PREFIX=${SAMPLE}/" \
	/FitHiChIP/configfile_BiasCorrection_CoverageBias \
	> /FitHiChIP/configfile_BiasCorrection_CoverageBias_${SAMPLE} && \
	bash /FitHiChIP/FitHiChIP_HiCPro.sh \
	-C /FitHiChIP/configfile_BiasCorrection_CoverageBias_${SAMPLE} \
	> ${FITHICHIPDIR}/${SAMPLE}.log 2>&1

done