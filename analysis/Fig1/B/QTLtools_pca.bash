#!/bin/bash

docker run --rm -v /:/wakasugi naotokubota/qtltools:1.3.1 bash -c "\
QTLtools pca \
--vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.vcf.gz \
--scale --center --maf 0.01 --distance 50000 \
--out /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/pca/merge.MAF001 && \
QTLtools pca \
--vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SVs_1KGP_pgGTs.438.MAF001.vcf.gz \
--scale --center --maf 0.01 --distance 50000 \
--out /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/pca/SVs_1KGP_pgGTs.438.MAF001 && \
QTLtools pca \
--vcf /wakasugi/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.sort.vcf.gz \
--scale --center --maf 0.01 --distance 50000 \
--out /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/pca/SV.SNP.sort\
"
