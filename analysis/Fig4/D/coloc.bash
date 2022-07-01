#!/bin/bash

####################### Must run in docker container! #########################
# docker run -v /:/wakasugi --rm -it naotokubota/coloc-locuscomparer:1.0 bash
###############################################################################

DATAPATH=/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt

# Make statistics file of puQTL for prmtr.64998
ID=prmtr.64998 && \
GENE=MSTRG.19786 && \
echo -e 'rs_id\tprmtr_id\tgene_id\tpval_nominal\tslope\tslope_se' \
> ${DATAPATH}/coloc/coloc_puQTL_${ID}.txt && \
zcat ${DATAPATH}/nominal_PEER25_full.txt.gz | \
awk -F" " -v OFS="\t" -v PRMTR=${ID} -v GENE=${GENE} '$1 == PRMTR{print $8,$1,GENE,$12,$14,$15}' \
>> ${DATAPATH}/coloc/coloc_puQTL_${ID}.txt

######################### Ulcertive colitis (GCST003043) ##########################

# Set Variables
GWAS=GCST003043
PROPOTION=0.37 # 12882 (Case) / 34652 (All)

# Make statistics file of GWAS for Ulcertive colitis (GCST003043)
echo -e 'rs_id\tpval\tslope\tslope_se\tgwas' \
> ${DATAPATH}/coloc/coloc_GWAS_${GWAS}.txt && \
zcat /wakasugi/mnt/data5/naoto/GWASCatalog/summary/*${GWAS}*.h.tsv.gz | \
awk -F"\t" -v OFS="\t" -v GWAS=${GWAS} '{print $2,$14,$21,$23,GWAS}' | \
sed -e '1d' \
>> ${DATAPATH}/coloc/coloc_GWAS_${GWAS}.txt
# Run coloc
Rscript coloc_vs_GWAS_cc.R \
${DATAPATH}/coloc/coloc_puQTL_${ID}.txt \
${DATAPATH}/coloc/coloc_GWAS_${GWAS}.txt \
${PROPOTION} \
${DATAPATH}/coloc/coloc_result_${ID}_${GWAS}.txt && \
echo -e "Done ${ID}-${GWAS}!"

####################################################################################