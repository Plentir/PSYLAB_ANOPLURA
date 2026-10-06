#!/bin/bash
#
# This script was originally authored by Daren Card (https://darencard.net/blog/2020-07-23-augustus-optimization/).
# For this study, we adapted it into the current wrapper script.
#
catalog="insecta_odb10"
species_name=<ROOT_DIR>
genome_name=<GENOME_FILE>
busco_result_dir=${species_name}/run_${catalog}

counter=1
cat ${catalog}.info.ranked.txt | grep -v "^#" | cut -f 1 | while read id;
do
    grep -w "${id}" complete_singlecopy.seqs.rename.cdhit.fasta
done | tr -d ">" | head -1334 | while read busco;
do
    cat ${busco_result_dir}/augustus_output/predicted_genes_rerun/${busco}.out.1 | sed "s/g1/g${counter}/g"
    ((counter=counter+1))
done | grep -v "^#" | \
    awk -F "\t" -v OFS="\t" '{split($1,a,/:/); print a[1], a[2], $2, $3, $4, $5, $6, $7, $8, $9 }' | \
    awk -F "\t" -v OFS="\t" '{split($2,a,/-/); print $1, $3, $4, a[1]+$5, a[1]+$6, $7, $8, $9, $10 }' | \
    sort -k1,1 -k4,4n > ${catalog}.info.ranked.complete.nonredundant.gff

$TOOLBOX/Augustus-3.5.0/scripts/gff2gbSmallDNA.pl ${catalog}.info.ranked.complete.nonredundant.gff \
    ${genome_name} 1000 \
    ${catalog}.info.ranked.complete.nonredundant.gb

$TOOLBOX/Augustus-3.5.0/scripts/randomSplit.pl ${catalog}.info.ranked.complete.nonredundant.gb 25

$TOOLBOX/Augustus-3.5.0/scripts/new_species.pl --species=${species_name}

# Initial run (Untrained state)
$TOOLBOX/Augustus-3.5.0/bin/etraining --species=${species_name} ${catalog}.info.ranked.complete.nonredundant.gb.train

$TOOLBOX/Augustus-3.5.0/bin/augustus --species=${species_name} ${catalog}.info.ranked.complete.nonredundant.gb.test | tee augustus_untrained.stat

# Optimization, then Re-run (Trained state)
$TOOLBOX/Augustus-3.5.0/scripts/optimize_augustus.pl --aug_exec_dir=$TOOLBOX/Augustus-3.5.0/bin --cpus=12 --kfold=24 --species=${species_name} ${catalog}.info.ranked.complete.nonredundant.gb.train

$TOOLBOX/Augustus-3.5.0/bin/etraining --species=${species_name} ${catalog}.info.ranked.complete.nonredundant.gb.train

$TOOLBOX/Augustus-3.5.0/bin/augustus --species=${species_name} ${catalog}.info.ranked.complete.nonredundant.gb.test | tee augustus_trained.stat