#!/bin/bash
#
# This script was originally authored by Daren Card (https://darencard.net/blog/2020-07-23-augustus-optimization/).
# For this study, we adapted it into the current wrapper script.
#
cat 
    <(echo "#group\tgene_name\torthodb_url\tevolutionary_rate\tmedian_exon_count\tstdev_exon_count\tmedian_protein_length\tstdev_protein_length")
    <(cat links_to_ODB10_insecta.txt | cut -f 1 | while read id;
        do 
            curl -s "https://v101.orthodb.org/group?id=${id}" | jq -r '. | [.data.public_id, .data.name, .url, .data.evolutionary_rate, .data.gene_architecture.exon_median_counts, .data.gene_architecture.exon_stdev_counts, .data.gene_architecture.protein_median_length, .data.gene_architecture.protein_stdev_length] | @tsv'
            sleep 1s
        done) > insecta_odb10.info.txt

cat <(cat insecta_odb10.info.txt | head -1) <(cat insecta_odb10.info.txt | grep -v "^#" | sort -t $'\t' -k5,5nr -k7,7nr -k4,4n) > insecta_odb10.info.ranked.txt
