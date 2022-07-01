#!/bin/bash

if [ ! -f /home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM_forpeer.tsv ]; then
	zcat /home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed.gz | \
	cut -f 7- | sed -e '1d' \
	> /home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM_forpeer.tsv
fi

docker run --rm -v /:/wakasugi broadinstitute/gtex_eqtl:V8 bash -c "\
Rscript /wakasugi/home/naoto/TSSchoiceQTL/scripts/peer.R \
/wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM_forpeer.tsv \
/wakasugi/home/naoto/TSSchoiceQTL/data/eQTL/txt/TMM_peerfactors.txt\
"
