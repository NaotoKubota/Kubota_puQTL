#!/bin/bash

### make TSS bed file ###

awk -F"\t" \
'$3 == "transcript"{print $9}' \
/home/naoto/data6/ENCODE/GM12878/CAGE/gencode.v24.primary_assembly.annotation.gtf | \
cut -d";" -f 1,2 | cut -d" " -f 2,4 | \
sed -e 's/"//g' -e 's/;//g' | awk -F" " -v OFS="\t" '{print $2,$1,$3}' \
> /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol.tsv

awk -F"\t" -v OFS="\t" \
'$3 == "transcript" && $7 == "+"{print $1,$4-1,$4,$7}\
$3 == "transcript" && $7 == "-"{print $1,$5-1,$5,$7}' \
/home/naoto/data6/ENCODE/GM12878/CAGE/gencode.v24.primary_assembly.annotation.gtf \
> /home/naoto/data6/ENCODE/GM12878/CAGE/ENST.bed

paste \
/home/naoto/data6/ENCODE/GM12878/CAGE/ENST.bed \
/home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol.tsv | \
sort -k1,1 -k2,2n \
> /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS.bed

rm -rf /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol.tsv
rm -rf /home/naoto/data6/ENCODE/GM12878/CAGE/ENST.bed
head /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS.bed

### cluster TSS (within 50 bp) and give IDs ###
cut -f 6 /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS.bed | awk '!a[$0]++' | while read line;
do

awk -F"\t" -v ENSG=${line} -v OFS="\t" '$6 == ENSG{print}' /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS.bed | \
bedtools cluster -d 50 -i stdin | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$4,$5,$6"_TSS_"$7}' \
>> /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS_cluster_tmp.bed

done

sort -k1,1 -k2,2n /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS_cluster_tmp.bed \
> /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS_cluster.bed

rm -rf /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS_cluster_tmp.bed

# fastp
fastp \
-i ENCFF000EWJ.fastq.gz \
-I ENCFF000EWX.fastq.gz \
-3 \
-o ENCFF000EWJ_cleaned.fastq.gz \
-O ENCFF000EWX_cleaned.fastq.gz \
-h fastp_report.html \
-q 30 -w 24

# Random sampling (50M reads)
seqkit sample -n 50000000 -s 20 ENCFF000EWJ_cleaned.fastq.gz > ENCFF000EWJ_cleaned_50M.fastq.gz
seqkit sample -n 50000000 -s 20 ENCFF000EWX_cleaned.fastq.gz > ENCFF000EWX_cleaned_50M.fastq.gz

# kallisto
## make index ##
kallisto index \
-i /home/naoto/data6/kallisto/homo_sapiens/gencode.v24.transcripts.idx \
/home/naoto/data6/kallisto/homo_sapiens/gencode.v24.transcripts.fa.gz

## quant ##
kallisto quant \
-t 12 \
-i /home/naoto/data6/kallisto/homo_sapiens/gencode.v24.transcripts.idx \
-o ENCFF000EWJ_ENCFF000EWX \
ENCFF000EWJ_cleaned_50M.fastq.gz \
ENCFF000EWX_cleaned_50M.fastq.gz

%%bash

# Sum TPM for each TSS
rm -rf /home/naoto/data6/ENCODE/GM12878/RNAseq/ENCFF000EWJ_ENCFF000EWX/TPM_TSS_cluster.tsv

awk -F"\t" -v OFS="\t" '{split($1,ID,"|"); print $5,ID[1]}' \
/home/naoto/data6/ENCODE/GM12878/RNAseq/ENCFF000EWJ_ENCFF000EWX/abundance.tsv | \
sed -s '1d' | awk -F"\t" -v OFS="\t" '{print $2,$1}' | sort | \
join -t $'\t' -1 5 -2 1 \
<(sort -k5 /home/naoto/data6/ENCODE/GM12878/CAGE/ENST_ENSG_Symbol_TSS_cluster.bed) \
- | cut -f 1,6,7 | \
awk -F"\t" -v OFS="\t" '{ if($0 != "") { a[$2]+=$3; } }END{for(i in a)print i,a[i];}' \
>> /home/naoto/data6/ENCODE/GM12878/RNAseq/ENCFF000EWJ_ENCFF000EWX/TPM_TSS_cluster.tsv
