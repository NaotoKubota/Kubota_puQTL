#!/bin/bash

awk -F"\t" -v OFS="\t" '$7 == "False" || $1 == "seqnames"' /home/naoto/TSSchoiceQTL/data/proActiv/proActiv_GEUVADIS_merged.bed | \
awk -F"\t" -v OFS="\t" '$1 != "chrX"' | \
awk -F"\t" -v OFS="\t" '$1 != "chrY"' | \
cut -f 1-6,9- | sed -e 's/_SJ.out//g' -e 's/seqnames/#Chr/g' | \
awk -F"\t" -v OFS="\t" '{split($4,id,"."); $4="prmtr."id[1]; print}' | \
sort -k1,1 -k2,2n \
> /home/naoto/TSSchoiceQTL/data/QTLtools/bed/proActiv_GEUVADIS_merged_QTLtools.bed

bgzip /home/naoto/TSSchoiceQTL/data/QTLtools/bed/proActiv_GEUVADIS_merged_QTLtools.bed && \
tabix -p bed /home/naoto/TSSchoiceQTL/data/QTLtools/bed/proActiv_GEUVADIS_merged_QTLtools.bed.gz

# chr22
awk -F"\t" -v OFS="\t" '$7 == "False" || $1 == "seqnames"' /home/naoto/TSSchoiceQTL/data/proActiv/proActiv_GEUVADIS_merged.bed | \
awk -F"\t" -v OFS="\t" '$1 != "chrX"' | \
awk -F"\t" -v OFS="\t" '$1 != "chrY"' | \
awk -F"\t" -v OFS="\t" '$1 == "seqnames" || $1 == "chr22"' | \
cut -f 1-6,9- | sed -e 's/_SJ.out//g' -e 's/seqnames/#Chr/g' | \
awk -F"\t" -v OFS="\t" '{split($4,id,"."); $4="prmtr."id[1]; print}' | \
sort -k1,1 -k2,2n \
> /home/naoto/TSSchoiceQTL/data/QTLtools/bed/proActiv_GEUVADIS_merged_QTLtools_chr22.bed

bgzip /home/naoto/TSSchoiceQTL/data/QTLtools/bed/proActiv_GEUVADIS_merged_QTLtools_chr22.bed && \
tabix -p bed /home/naoto/TSSchoiceQTL/data/QTLtools/bed/proActiv_GEUVADIS_merged_QTLtools_chr22.bed.gz

