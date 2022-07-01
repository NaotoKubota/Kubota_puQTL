#!/bin/bash

############################# Must run in docker container! #################################
# docker run --rm -it -v /:/wakasugi quay.io/biocontainers/plink:1.90b6.21--h779adbc_1 bash
#############################################################################################

zcat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz | \
awk '$23 == 1{print $8,$2}' | awk '!a[$0]++' \
> /home/naoto/TSSchoiceQTL/data/puQTL_PEER25_besthit_chr_rs.txt

# calculate LD
rm -rf /home/naoto/TSSchoiceQTL/data/ld_puQTL_besthit_PEER25.txt && \
cat /home/naoto/TSSchoiceQTL/data/puQTL_PEER25_besthit_chr_rs.txt | \
while read line
do

    ID=`echo ${line} | cut -d' ' -f 1` && \
    CHR=`echo ${line} | cut -d' ' -f 2` && \
    echo -e "START ${ID}" && \
    plink \
    --vcf /rhome/naotok/bigdata/GEUVADIS/vcf/SV.SNP.${CHR}.sort.vcf.gz \
    --r2 \
    --ld-snp ${ID} \
    --ld-window-kb 1000 \
    --ld-window 99999 \
    --ld-window-r2 0.5 \
    --threads 50 && \
    cat plink.ld | \
    sed -e 's/  */\t/g' -e 's/^\t//g' | \
    cut -f 3,6,7 | \
    sed -e 1d \
    >> /home/naoto/TSSchoiceQTL/data/ld_puQTL_besthit_PEER25.txt && \
    rm -rf plink.ld plink.nosex plink.log && \
    echo -e "Done ${ID}"

done