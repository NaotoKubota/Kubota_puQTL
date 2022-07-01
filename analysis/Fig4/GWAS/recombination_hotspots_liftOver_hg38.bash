#!/bin/bash

# concat
cat /mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh37_chr*.txt | \
grep -v "Chromosome" | \
grep -v "_" | \
awk -F"\t" -v OFS="\t" '{print $1,$2-1000,$2+1000,$3}' | \
sort -k1,1 -k2,2n \
> /mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh37.bed

# liftOver (hg19 to hg38)
/home/naoto/src/liftOver \
/mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh37.bed \
/mnt/data5/naoto/GEUVADIS/vcf/hg19ToHg38.over.chain.gz \
/mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh38.bed \
/mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh37_notConverted.bed

# cM > 10 and merge
cat /mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh38.bed | \
awk -F"\t" '$4 > 10' | \
sort -k1,1 -k2,2n | \
bedtools merge -i stdin \
> /mnt/data5/naoto/GEUVADIS/vcf/genetic_map_GRCh38_cM10_merge.bed
