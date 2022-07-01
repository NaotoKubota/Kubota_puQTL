#!bin/bash

zcat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz | \
cut -d " " -f 1,2 | awk '!a[$0]++' \
> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene.txt && \
split -l 200 -d /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene.txt \
/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene- && \
rm -rf /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene.txt


function forsusie () {

	cat ${1} | while read line
	do

		PHENOTYPE=`echo ${line} | cut -d" " -f 1`
		CHR=`echo ${line} | cut -d" " -f 2`

		echo ${PHENOTYPE}
		echo ${CHR}

		# Get summary from nominal pass output
		zcat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/nominal_PEER25_full.txt.gz | 
		awk -v ID=${PHENOTYPE} -v OFS="\t" '$1 == ID{print $8,$14,$15}' | \
		awk '!a[$1]++' \
		> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${PHENOTYPE}.txt && \
		cat /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${PHENOTYPE}.txt | \
		cut -f 1 \
		> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${PHENOTYPE}.variant.txt

		# Get genotype
		bcftools query \
		--include ID=@/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${PHENOTYPE}.variant.txt \
		-f '%ID[\t%GT]\n' \
		/mnt/data5/naoto/GEUVADIS/vcf/SV.SNP.${CHR}.sort.vcf.gz | \
		sed \
		-e 's@0/0@0@g' \
		-e 's@0|0@0@g' \
		-e 's@0/1@1@g' \
		-e 's@0|1@1@g' \
		-e 's@1/0@1@g' \
		-e 's@1|0@1@g' \
		-e 's@1/1@2@g' \
		-e 's@1|1@2@g' | \
		awk '!a[$1]++' \
		> /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/${PHENOTYPE}.genotype.txt

	done


}


# Parallel run
for file in `ls /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene-*`;
do
	
	echo ${file}

	forsusie ${file} &

done

wait

rm -rf /home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_gene-*

echo "Done!"
