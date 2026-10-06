#!/bin/bash
#
# This script was originally authored by Daren Card (https://darencard.net/blog/2020-07-23-augustus-optimization/).
# For this study, we adapted it into the current wrapper script.
#
catalog="insecta_odb10"
busco_result_dir=<ROOT_DIR>/run_${catalog}

rm -f complete_singlecopy.aafiles.txt complete_singlecopy.seqs.rename.fasta complete_singlecopy.seqs.rename.fasta.cdhit
processed_job=1

grep -v "^#" ${catalog}.info.ranked.txt | cut -f 1 | while read id;
do
    echo "Job $processed_job is now running."
    file=${busco_result_dir}/augustus_output/extracted_proteins/"${id}".faa.1

	status=`awk -v ID="$id" '{ if ($1 == ID) print $2 }' ${busco_result_dir}/full_table.tsv | sort | uniq`
	count=`grep -c "^>" $file`

    echo "${id}: ${status} / ${count}"
    
    if [ "$count" = "1" ] && [ "$status" = "Complete" ]; then
        echo "==> BUSCO ID $id has been retained."
        cat $file | $TOOLBOX/seqkit fx2tab | awk -v ID="$id" -F "\t" '{ print ID "\t" $2 }' | $TOOLBOX/seqkit tab2fx >> complete_singlecopy.seqs.rename.fasta
    fi
    ((processed_job=processed_job+1))
done

$TOOLBOX/cdhit-4.8.1/cd-hit -o complete_singlecopy.seqs.rename.fasta.cdhit -c 0.8 -i complete_singlecopy.seqs.rename.fasta -p 1 -d 0 -b 3 -T 0 -M 16000
mv complete_singlecopy.seqs.rename.fasta.cdhit complete_singlecopy.seqs.rename.cdhit.fasta
