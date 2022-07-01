#!/bin/bash

####################### Must run in docker container! #########################
# docker run -v /:/wakasugi --rm -it naotokubota/coloc-locuscomparer:1.0 bash
###############################################################################

DATAPATH=/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt

# Make statistics file of puQTL for prmtr.109751
ID=prmtr.109751 && \
GENE=MSTRG.33276 && \
echo -e 'rs_id\tprmtr_id\tgene_id\tpval_nominal\tslope\tslope_se' \
> ${DATAPATH}/coloc/coloc_puQTL_${ID}.txt && \
zcat ${DATAPATH}/nominal_PEER25_full.txt.gz | \
awk -F" " -v OFS="\t" -v PRMTR=${ID} -v GENE=${GENE} '$1 == PRMTR{print $8,$1,GENE,$12,$14,$15}' \
>> ${DATAPATH}/coloc/coloc_puQTL_${ID}.txt


######################## Diverticular disease (GCST006479) ########################

# Set variables
GWAS=GCST006479
PROPOTION=0.067 # 27444 (Case) / 409728 (All)

# Make statistics file of GWAS for Diverticular disease (GCST006479)
echo -e 'rs_id\tpval\tslope\tslope_se\tgwas' \
> ${DATAPATH}/coloc/coloc_GWAS_${GWAS}.txt && \
zcat /wakasugi/mnt/data5/naoto/GWASCatalog/summary/clinical_c_K57/imputed.allWhites.clinical_c_K57.chr9.csv.gz | \
awk -F" " -v OFS="\t" -v GWAS=${GWAS} '{print $1,$6,$4,$5,GWAS}' | \
sed -e '1d' \
>> ${DATAPATH}/coloc/coloc_GWAS_${GWAS}.txt
# Run coloc
Rscript coloc_vs_GWAS_cc.R \
${DATAPATH}/coloc/coloc_puQTL_${ID}.txt \
${DATAPATH}/coloc/coloc_GWAS_${GWAS}.txt \
${PROPOTION} \
${DATAPATH}/coloc/coloc_result_${ID}_${GWAS}.txt && \
echo -e "Done ${ID}-${GWAS}!"

###################################################################################