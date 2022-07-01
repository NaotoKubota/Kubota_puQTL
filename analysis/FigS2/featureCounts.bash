#!/bin/bash

DATAPATH=/mnt/data4/naoto/GEUVADIS/STAR

# featureCounts
featureCounts -p -B -t exon -g gene_id \
-a /mnt/data4/naoto/GEUVADIS/STAR/mapping/gtf/GEUVADIS_merge.strands.gtf \
-o ${DATAPATH}/featureCounts/all_counts_stringtie.txt \
-T 16 \
${DATAPATH}/mapping/*/2nd/*.sort.bam

# filter
python expression_filter.py

# TMM
docker run --rm -v /:/wakasugi broadinstitute/gtex_eqtl:V8 bash -c "\
Rscript /wakasugi/home/naoto/TSSchoiceQTL/scripts/eQTL/TMM.R\
" && \
python make_input_table.py && \
bgzip /home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed && \
tabix -p bed /home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed.gz
