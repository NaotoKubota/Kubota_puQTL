#!/bin/bash

# cp /home/naoto/data4/GEUVADIS/STAR/mapping/*/2nd/*_SJ.out.tab /home/naoto/data4/GEUVADIS/STAR/junctions/
# Remove empty files
find /home/naoto/data4/GEUVADIS/STAR/junctions/ -empty -delete

# drop transcripts without strand information
grep -v -P '\t\.\t\.\t' /mnt/data4/naoto/GEUVADIS/STAR/mapping/gtf/GEUVADIS_merge.gtf \
> /mnt/data4/naoto/GEUVADIS/STAR/mapping/gtf/GEUVADIS_merge.strands.gtf

# Estimating promoter activity by proActiv
docker run --rm -v /:/wakasugi naotokubota/proactiv:1.1.18 bash -c "\
Rscript /wakasugi/home/naoto/TSSchoiceQTL/scripts/proActiv/proActiv_GEUVADIS.R" && \
cat \
<(echo -e "seqnames\tstart\tend\tpromoterId\tgeneId\tstrand\tinternalPromoter\tpromoterPosition\ttxId") \
<(cut -f 2- /home/naoto/TSSchoiceQTL/data/proActiv/result_rowData_GEUVADIS.tsv | sed -e 1d | awk -F"\t" -v OFS="\t" '{print $3,$4,$4+1,$1,$2,$5,$6,$7,$8}') | \
awk -F"\t" -v OFS="\t" '$1 ~ "chr" || $1 == "seqnames"' | \
paste \
- \
<(cat <(head -n 1 /home/naoto/TSSchoiceQTL/data/proActiv/result_absolutePromoterActivity_GEUVADIS.tsv) <(cut -f 2- /home/naoto/TSSchoiceQTL/data/proActiv/result_absolutePromoterActivity_GEUVADIS.tsv | sed -e 1d)) \
> /home/naoto/TSSchoiceQTL/data/proActiv/proActiv_GEUVADIS_merged_stringtie.bed && \
sort -k1,1 -k2,2n /home/naoto/TSSchoiceQTL/data/proActiv/proActiv_GEUVADIS_merged_stringtie.bed \
> /home/naoto/TSSchoiceQTL/data/proActiv/proActiv_GEUVADIS_merged_sorted_stringtie.bed

# Drop internal promoters
awk -F"\t" -v OFS="\t" '$7 == "FALSE" || $1 == "seqnames"' /home/naoto/TSSchoiceQTL/data/proActiv/proActiv_GEUVADIS_merged_stringtie.bed | \
cut -f 1-6,10- | sed -e 's/_SJ.out//g' -e 's/seqnames/#Chr/g' | \
awk -F"\t" -v OFS="\t" '{$4="prmtr."$4; print}' | \
sort -k1,1 -k2,2n \
> /home/naoto/TSSchoiceQTL/data/puQTL/bed/proActiv_GEUVADIS_merged_stringtie_QTLtools.bed && \
bgzip /home/naoto/TSSchoiceQTL/data/puQTL/bed/proActiv_GEUVADIS_merged_stringtie_QTLtools.bed && \
tabix -p bed /home/naoto/TSSchoiceQTL/data/puQTL/bed/proActiv_GEUVADIS_merged_stringtie_QTLtools.bed.gz

# Remove the promoters of which the 4th quantile of score is 0
python proActiv_filter.py && \
bgzip /home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed && \
tabix -p bed /home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed.gz
