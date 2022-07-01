#!/bin/bash

# rmdup
# SNP
bcftools view \
--threads 20 \
-m2 -M2 \
--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
/mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.vcf.gz | \
bcftools norm \
--threads 20 \
--rm-dup all | \
bgzip -c \
> /mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.biallelic.rmdup.vcf.gz && \
tabix -p vcf /mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.biallelic.rmdup.vcf.gz
# SV
bcftools view \
--threads 20 \
-m2 -M2 \
--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
/mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.vcf.gz | \
bcftools norm \
--threads 20 \
--rm-dup all | \
bgzip -c \
> /mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.biallelic.rmdup.vcf.gz && \
tabix -p vcf /mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.biallelic.rmdup.vcf.gz

# concat SNP & SV for each chromosome
for i in chr{1..22} chrX chrY
do

	echo ${i} && \
	bcftools view \
	--threads 20 \
	/mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.biallelic.rmdup.vcf.gz \
	--regions ${i} | \
	bgzip -c \
	> /mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.biallelic.rmdup.${i}.vcf.gz && \
	bcftools view \
	--threads 20 \
	/mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.biallelic.rmdup.vcf.gz \
	--regions ${i} | \
	bgzip -c \
	> /mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.biallelic.rmdup.${i}.vcf.gz && \
	bcftools concat \
	--threads 20 \
	/mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.biallelic.rmdup.${i}.vcf.gz \
	/mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.biallelic.rmdup.${i}.vcf.gz | \
	bgzip -c \
	> /mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.${i}.vcf.gz

done

# sort
for i in chr{1..22} chrX chrY
do

	bcftools sort /mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.${i}.vcf.gz \
	--temp-dir /mnt/data5/naoto/GEUVADIS/vcf/ | \
	bgzip -c \
	> /mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.${i}.sort.vcf.gz

done

# concat all chromosomes
bcftools concat \
--threads 20 \
/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.*.sort.vcf.gz | \
bcftools sort \
--temp-dir /mnt/data5/naoto/GEUVADIS/vcf/ | \
bgzip -c \
> /mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz
