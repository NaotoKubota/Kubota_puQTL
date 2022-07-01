#!/bin/bash

############################## Must run in docker container ##############################
# docker run --rm -it -v /:/wakasugi naotokubota/coloc-locuscomparer:1.0 bash
##########################################################################################

zcat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz | \
cut -d " " -f 1,2 | awk '!a[$0]++' | sort \
> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene.txt && \
split -l 200 -d /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene.txt \
/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene- && \
rm -rf /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene.txt

# All phenoype
zcat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz | \
cut -d " " -f 1 | awk '!a[$0]++' | sort \
> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_all.txt && \
# Processed phenoype
ls /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/*result.txt | \
cut -d"." -f -2 | cut -d"/" -f 9 | sort \
> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_processed.txt && \
# Unprocessed phenoype
comm -23 \
/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_all.txt \
/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_processed.txt \
> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_unprocessed.txt && \
join \
/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene.txt \
/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_unprocessed.txt \
> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_unprocessed_ID_chr.txt

# Remove duplicated variants from genotype and phenotype table
cat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_unprocessed.txt | while read line
do

    # genotype
    cat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.genotype.txt | \
    awk '!a[$1]++' \
    > /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.genotype.2.txt && \
    rm -rf /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.genotype.txt && \
    mv \
    /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.genotype.2.txt \
    /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.genotype.txt
    # phenotype
    cat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.txt | \
    awk '!a[$1]++' \
    > /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.2.txt && \
    rm -rf /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.txt && \
    mv \
    /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.2.txt \
    /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${line}.txt

done

# SuSiE
function SuSiE () {

    cat ${1} | while read line
    do

        echo -e "START ${line}"
        Rscript susie.R ${line}

    done

}

# Parallel run
for file in `ls /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene-*`;
do
	
	echo ${file}

	SuSiE ${file} &

done

wait


# Single run for unprocessed phenotypes
SuSiE /wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/susie_unprocessed_ID_chr.txt
