#!/bin/bash

# HiC-Pro
docker run -v /:/wakasugi --rm nservant/hicpro bash -c "\
/HiC-Pro_3.0.0/bin/utils/digest_genome.py \
-r ^GATC -o /HiC-Pro_3.0.0/annotation/DpnII_resfrag_hg38.bed \
/wakasugi/mnt/data6/naoto/twins1/RefGenome/hg38.fa && \
cp \
/wakasugi/mnt/data6/naoto/twins1/RefGenome/hg38_genomeSize.txt \
/HiC-Pro_3.0.0/annotation/chrom_hg38.sizes && \
sed \
-e 's/N_CPU\ =\ 2/N_CPU\ =\ 24/' \
-e 's/SORT_RAM\ =\ 768M/SORT_RAM\ =\ 1000M/' \
-e 's/BOWTIE2_IDX_PATH\ =/BOWTIE2_IDX_PATH\ =\ \/wakasugi\/mnt\/data6\/common\/Mapping_Index\/bowtie2_hg38\//' \
-e 's/REFERENCE_GENOME\ =\ hg19/REFERENCE_GENOME\ =\ hg38/' \
-e 's/GENOME_SIZE\ =\ chrom_hg19\.sizes/GENOME_SIZE\ =\ \/HiC-Pro_3\.0\.0\/annotation\/chrom_hg38\.sizes/' \
-e 's/GENOME_FRAGMENT\ =\ HindIII_resfrag_hg19\.bed/GENOME_FRAGMENT\ =\ \/HiC-Pro_3\.0\.0\/annotation\/DpnII_resfrag_hg38\.bed/' \
-e 's/LIGATION_SITE\ =\ AAGCTAGCTT/LIGATION_SITE\ =\ GATCGATC/' \
/HiC-Pro_3.0.0/config-hicpro.txt \
> /wakasugi/home/naoto/TSSchoiceQTL/data/HiChIP/config-hicpro-DpnII-hg38.txt && \
/HiC-Pro_3.0.0/bin/HiC-Pro \
-i /wakasugi/mnt/data4/naoto/puQTL/HiChIP/fastq/ \
-o /wakasugi/mnt/data4/naoto/puQTL/HiChIP/HiC-Pro_results_hg38 \
-c /wakasugi/home/naoto/TSSchoiceQTL/data/HiChIP/config-hicpro-DpnII-hg38.txt && \
BAMPATH=/wakasugi/mnt/data4/naoto/puQTL/HiChIP/HiC-Pro_results_hg38/bowtie_results/bwt2 && \
for SAMPLE in SRR5831489 SRR5831490; do \
samtools \
sort \
-@ 24 \
-o ${BAMPATH}/${SAMPLE}/${SAMPLE}_filtered_hg38.bwt2pairs.sort.bam \
${BAMPATH}/${SAMPLE}/${SAMPLE}_filtered_hg38.bwt2pairs.bam && \
samtools index -@ 24 ${BAMPATH}/${SAMPLE}/${SAMPLE}_filtered_hg38.bwt2pairs.sort.bam; \
done\
"