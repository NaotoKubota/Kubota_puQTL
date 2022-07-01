#!/bin/bash

ls /home/naoto/data4/GEUVADIS/fastq/fastp_json/ | \
grep -e HG -e NA | \
sed -e 's/_fastp\.json//g' \
> /home/naoto/TSSchoiceQTL/data/txt/438_ID.txt && \
bcftools view \
-S /home/naoto/TSSchoiceQTL/data/txt/438_ID.txt \
/mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.vcf.gz | \
grep -v NO_VALID_GT | \
bgzip -c \
> /mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.vcf.gz && \
bcftools +fill-tags \
/mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.vcf.gz \
-- -t AF | \
bcftools view --exclude 'AF[0]<0.01 | AF[0]>0.99' | \
bgzip -c \
> /mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.vcf.gz && \
tabix -p vcf /mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.vcf.gz
