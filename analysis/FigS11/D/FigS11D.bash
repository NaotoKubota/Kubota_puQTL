#!/bin/bash

############################# Must run in docker container! #################################
# docker run --rm -it -v /:/wakasugi quay.io/biocontainers/plink:1.90b6.21--h779adbc_1 bash
#############################################################################################

# rs10229453
plink \
--vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.chr7.sort.vcf.gz \
--r2 \
--ld-snp rs10229453 \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--threads 24 && \
cat plink.ld | \
sed -e 's/  */\t/g' -e 's/^\t//g' | \
cut -f 4-7 | \
sed -e 1d \
> /wakasugi/home/naoto/TSSchoiceQTL/data/txt/coloc/ld_rs10229453.txt && \
rm -rf plink.ld plink.nosex plink.log
