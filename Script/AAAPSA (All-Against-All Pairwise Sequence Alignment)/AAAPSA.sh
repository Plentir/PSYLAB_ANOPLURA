#!/bin/bash
seqkit_dir=<PATH_TO_SEQKIT>

help(){
    echo "Calculate pairwise sequence identity by all possible pairs"
    echo "Usage: ./AAAPSA.sh -i <set_of_sequences_1> -I <set_of_sequences_2> -m <dna|rna|protein> -g"
    echo "-i1, -i2: Input FASTA file. Be aware the inputs must follow the format below."
    echo "          >ORTHOGROUP_IDENTIFIER:SPECIES_NAME:ACCESSION_OR_YOUR_OWN_TITLE"
    echo "          ATTTACGGACGATTA"
    echo ""
    echo "          <Example>"
    echo "          >OG0000001:human:NP_116614.1"
    echo "          (sequence)"
    echo "          * Currently input2 option is not available."
    echo ""
    echo "-g : Use global alignment if the flag exists, otherwise use local alignment."
    echo "     input file 'i2' will be ignored when using global alignment."
    echo "-m : Molecular type of sequences."
    echo "-s : Choose scoring matrix. Available matrices can be found here."
}

no_args="true"
in1=false
in2=false
isglobal=false
while getopts "i:I:m:s:gh" opt; do
    case $opt in
        i)
            in1=$OPTARG;;
        I)
            in2=$OPTARG;;
        m)
            mol=$OPTARG;;
        g)
            isglobal=true;;
        s)
            matrix=$OPTARG;;
        h)
            help
            exit 0;;
        *)
            help
            exit 0;;
    esac
    no_args="false"
done

if [ $no_args == "true" ] ; then
    help
    exit 0
fi

if false ; then
mkdir -p AAAPSA_temp/"${in1}.split"

################################################################################ input 1
grep ">" ${in1} | cut -c 2- > AAAPSA_temp/${in1}.molecules.acc
cut -d: -f1 AAAPSA_temp/${in1}.molecules.acc | sort | uniq > AAAPSA_temp/${in1}.orthogroups.txt
cut -d: -f2 AAAPSA_temp/${in1}.molecules.acc | sort -f | uniq -i > AAAPSA_temp/${in1}.species.txt

end="$(cat AAAPSA_temp/${in1}.species.txt | wc -l)"
seq $end | sed -e "s/ /\n/g" > AAAPSA_temp/${in1}.species_ids.txt
paste AAAPSA_temp/${in1}.species.txt AAAPSA_temp/${in1}.species_ids.txt > AAAPSA_temp/${in1}.species_maps.txt

"$seqkit_dir"/seqkit replace -p ":(.*):" -r ":{kv}:" -k AAAPSA_temp/${in1}.species_maps.txt ${in1} | "$seqkit_dir"/seqkit sort -n -o AAAPSA_temp/${in1}.rename

while read og; do
    "$seqkit_dir"/seqkit grep -n -t "${mol}" -r -p "${og}:" AAAPSA_temp/${in1}.rename > AAAPSA_temp/${in1}.split/"${og}.fasta";
done < AAAPSA_temp/${in1}.orthogroups.txt

fi
"$seqkit_dir"/seqkit stats -Tab -o AAAPSA_temp/sequence_split_stats_${in1}.tsv AAAPSA_temp/${in1}.split/*.fasta

echo "OrthoGroup	Qry.Acc	Sbj.Acc	Qry.Spec	Sbj.Spec	Pct.ID	Pct.Sim	Pct.Gaps	nGaps	Qry.Start	Qry.End	Sbj.Start	Sbj.End	Qry.Cov	Sbj.Cov	Score	Program	GapOpen	GapExt.	Matrix" > result_${in1}.tsv

################################################################################ input 2
if [ -r ${in2} ] ; then
    mkdir -p AAAPSA_temp/"${in2}.split"

    grep ">" ${in2} | cut -c 2- > AAAPSA_temp/${in2}.molecules.acc
    cut -d: -f1 AAAPSA_temp/${in2}.molecules.acc | sort | uniq > AAAPSA_temp/${in2}.orthogroups.txt
    cut -d: -f2 AAAPSA_temp/${in2}.molecules.acc | sort -f | uniq -i > AAAPSA_temp/${in2}.species.txt
    
    end="$(cat AAAPSA_temp/${in2}.species.txt | wc -l)"
    seq $end | sed -e "s/ /\n/g" > AAAPSA_temp/${in2}.species_ids.txt
    paste AAAPSA_temp/${in2}.species.txt AAAPSA_temp/${in2}.species_ids.txt > AAAPSA_temp/${in2}.species_maps.txt
    
    "$seqkit_dir"/seqkit replace -p ":(.*):" -r ":{kv}:" -k AAAPSA_temp/${in2}.species_maps.txt ${in2} | "$seqkit_dir"/seqkit sort -n -o AAAPSA_temp/${in2}.rename
    
    while read og; do
        "$seqkit_dir"/seqkit grep -n -t "${mol}" -r -p "${og}:" AAAPSA_temp/${in2}.rename > AAAPSA_temp/${in2}.split/"${og}.fasta";
    done < AAAPSA_temp/${in2}.orthogroups.txt

    "$seqkit_dir"/seqkit stats -Tab -o AAAPSA_temp/sequence_split_stats_${in2}.tsv AAAPSA_temp/${in2}.split/*.fasta

    echo "OrthoGroup	Qry.Acc	Sbj.Acc	Qry.Spec	Sbj.Spec	Pct.ID	Pct.Sim	Pct.Gaps	nGaps	Qry.Start	Qry.End	Sbj.Start	Sbj.End	Qry.Cov	Sbj.Cov	Score	Program	GapOpen	GapExt.	Matrix" > result_${in1}_${in2}.tsv
fi

#################################################################################

rm AAAPSA_temp/*.species_ids.txt AAAPSA_temp/*.species.txt

if [ -r ${in2} ] ; then
    for f1 in AAAPSA_temp/"${in1}.split"/*.fasta; do
        for f2 in AAAPSA_temp/"${in2}.split"/*.fasta; do
            if [ "${in1}" -eq "${in2}" ] ; then
                echo "Now processing '${f1}' * '${f2}'";
                python AAAPSA_mp.py ${f1} ${f2} ${isglobal} ${mol} >> result_${in1}_${in2}.tsv;
            fi
        done;
    done;
else 
    for f in AAAPSA_temp/"${in1}.split"/*.fasta; do
        echo "Now processing '${f}'";
        python AAAPSA_mp.py ${f} ${f} ${isglobal} ${mol} >> result_${in1}.tsv;
    done;
fi

echo "Done."
