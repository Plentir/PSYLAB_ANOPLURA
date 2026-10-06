#!/bin/bash
#
# This script was originally authored by Daren Card (https://darencard.net/blog/2022-07-09-genome-repeat-annotation/).
# For this study, we adapted it into the current wrapper script.
# The 'repclassifier.sh' can be found in the article linked above.
#
scripts=<PATH_TO_repclassifier.sh>
threads=12
workdir=`pwd`
genome_name=<GENOME_FILE_NAME>
spec_name_6char=<ALIAS_OF_GENOME>

class="Neoptera"      # Variable 'class' is used for RepBase-based identification of unknown repeats
order="anoplura"      # Variable 'order' is used for repeatmasker


## Self logging of script
cat $0 > ${workdir}/RMPipe_`date +%y%m%d_%Hh%Mm%Ss`.cmd
echo "+ Variables"
echo "| Workdir=${workdir}"
echo "| Genome Name=${genome_name}"
echo "| Species Name (Abbr.)=${spec_name_6char}"
echo "+----------------------------------------"

## Preparation - Divide RepeatModeler 'families.fa' file into known (identified) and unknown repeats
echo "################### Now Preparing Data For RepeatMasker Rounds. ###################"
cat ${genome_name}-families.fa |
    ${TOOLBOX}/seqkit fx2tab |
    awk -v spec=${spec_name_6char} '{ print spec"_"$0 }' |
    ${TOOLBOX}/seqkit tab2fx > ${genome_name}-families.prefix.fa

cat ${genome_name}-families.prefix.fa |
    ${TOOLBOX}/seqkit fx2tab |
    grep -v "Unknown" |
    ${TOOLBOX}/seqkit tab2fx > ${genome_name}-families.prefix.fa.known

cat ${genome_name}-families.prefix.fa |
    ${TOOLBOX}/seqkit fx2tab |
    grep "Unknown" |
    ${TOOLBOX}/seqkit tab2fx > ${genome_name}-families.prefix.fa.unknown

# RepBase DB-based identification of unknown repeats (optional)
${scripts}/repclassifier.sh -t ${threads} -d ${class} \
    -u ${genome_name}-families.prefix.fa.unknown \
    -k ${genome_name}-families.prefix.fa.known \
    -a ${genome_name}-families.prefix.fa.known -o repbase_${class}_self

# Repeated identification (optional)
#${scripts}/repclassifier.sh -t 12 \
#    -u round-1_Repbase${class}-Self/round-1_Repbase${class}-Self.unknown \
#    -k round-1_Repbase${class}-Self/round-1_Repbase${class}-Self.known \
#    -a round-1_Repbase${class}-Self/round-1_Repbase${class}-Self.known -o round-2_Self

rm ${genome_name}-families.prefix.fa.known ${genome_name}-families.prefix.fa.unknown ${genome_name}-families.prefix.fa


## Round 1 - Find repeat sequences in de novo method
echo "################### Start RepeatMasker Round 1. ###################"
result_dir=r1.simple_repeats
${TOOLBOX}/RepeatMasker-4.1.6/RepeatMasker -pa ${threads} -s -a -e crossmatch -noint -xsmall -xm -gff \
    -species "${order}" \
    -dir ${workdir}/${result_dir} \
    ${workdir}/${genome_name}.fasta > ${workdir}/repmasker_round1.log

rename -E "s/fasta/simple_mask/" ${workdir}/${result_dir}/*
mv ${workdir}/${result_dir}/${genome_name}.simple_mask.masked ${workdir}/${result_dir}/${genome_name}.simple_mask.masked.fasta


## Round 2 - Find repeat sequences using repbase
echo "################### Start RepeatMasker Round 2. ###################"
input_dir=${result_dir}
result_dir=r2.${order}_repeats
${TOOLBOX}/RepeatMasker-4.1.6/RepeatMasker -pa ${threads} -s -a -e crossmatch -nolow -xsmall -xm -gff \
    -species ${order} \
    -dir ${workdir}/${result_dir} \
    ${workdir}/${input_dir}/${genome_name}.simple_mask.masked.fasta > ${workdir}/repmasker_round2.log

rename -E "s/simple_mask.masked.fasta/${order}_mask/" ${workdir}/${result_dir}/*
mv ${workdir}/${result_dir}/${genome_name}.${order}_mask.masked ${workdir}/${result_dir}/${genome_name}.${order}_mask.masked.fasta


## Round 3 - Find repeat sequences using identified repeats of Repeatmodeler predictions
echo "################### Start RepeatMasker Round 3. ###################"
input_dir=${result_dir}
result_dir=r3.modeler_known_repeats
${TOOLBOX}/RepeatMasker-4.1.6/RepeatMasker -pa ${threads} -s -a -e crossmatch -nolow -xsmall -xm -gff \
    -lib ${workdir}/repbase_${class}_self/repbase_${class}_self.known \
    -dir ${workdir}/${result_dir} \
    ${workdir}/${input_dir}/${genome_name}.${order}_mask.masked.fasta > ${workdir}/repmasker_round3.log

rename -E "s/${order}_mask.masked.fasta/modeler_known_mask/" ${workdir}/${result_dir}/*
mv ${workdir}/${result_dir}/${genome_name}.modeler_known_mask.masked ${workdir}/${result_dir}/${genome_name}.modeler_known_mask.masked.fasta


## Round 4 - Find repeat sequences using unknown repeats of Repeatmodeler predictions
echo "################### Start RepeatMasker Round 4. ###################"
input_dir=${result_dir}
result_dir=r4.modeler_unknown_repeats
${TOOLBOX}/RepeatMasker-4.1.6/RepeatMasker -pa ${threads} -s -a -e crossmatch -nolow -xsmall -xm -gff \
    -lib ${workdir}/repbase_${class}_self/repbase_${class}_self.unknown \
    -dir ${workdir}/${result_dir} \
    ${workdir}/${input_dir}/${genome_name}.modeler_known_mask.masked.fasta > ${workdir}/repmasker_round4.log

rename -E "s/modeler_known_mask.masked.fasta/modeler_unknown_mask/" ${workdir}/${result_dir}/*
mv ${workdir}/${result_dir}/${genome_name}.modeler_unknown_mask.masked ${workdir}/${result_dir}/${genome_name}.modeler_unknown_mask.masked.fasta


## Finish - Merge results of round 1 to 4
echo "################### Start Merging Results From RepeatMasker Round 1 to 4. ###################"
mkdir -p r5.all_merged
cat r[1-4].*/*.cat.gz > r5.all_merged/${genome_name}.merged_mask.cat.gz
cat r[1-4].*/*.align > r5.all_merged/${genome_name}.merged_mask.align

cat r1.simple_repeats/${genome_name}.simple_mask.out \
	<(tail -n +4 r2.${order}_repeats/${genome_name}.${order}_mask.out) \
	<(tail -n +4 r3.modeler_known_repeats/${genome_name}.modeler_known_mask.out) \
	<(tail -n +4 r4.modeler_unknown_repeats/${genome_name}.modeler_unknown_mask.out) \
	> r5.all_merged/${genome_name}.merged_mask.out

${TOOLBOX}/RepeatMasker-4.1.6/ProcessRepeats -a -species ${order} -xsmall -xm -gff \
    r5.all_merged/${genome_name}.merged_mask.cat.gz > r5.all_merged/${genome_name}.merged_mask.cat.gz.log


## Finish - Generate softmasked genome
bedtools maskfasta -soft -fi ${genome_name}.fasta -bed r5.all_merged/${genome_name}.merged_mask.out.gff \
    -fo r5.all_merged/${genome_name}.merged_mask.soft.fasta

echo "################### RepeatMasker Pipeline Has Been Done. ###################"
