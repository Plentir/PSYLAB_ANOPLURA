if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <paired_read1> <paired_read2>"
    exit 1
fi


toolbox_path=<PATH_TO_BBTOOLS_AND_LIGHTER>
accno=`echo $1 | sed -E "s/([^\.^_]+).+/\1/"`


echo "########################################################################################" > pipeline_$accno.log
echo "# Preprocessing Pipeline Started                                                       #" >> pipeline_$accno.log
echo "########################################################################################" >> pipeline_$accno.log
echo "" >> pipeline_$accno.log


echo "" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "       + BBTools: Clumpify (Group Read Mates and Remove Duplication) Started    +" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "" >> pipeline_$accno.log

$toolbox_path/bbmap/clumpify.sh in=$1 in2=$2 out=${accno}_1.clumped.fastq out2=${accno}_2.clumped.fastq groups=auto pigz dedupe 1>>pipeline_$accno.log 2>&1


echo "" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "       + BBTools: BBDuk (Remove Contamination of Reads) Started                 +" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "" >> pipeline_$accno.log

$toolbox_path/bbmap/bbduk.sh in=${accno}_1.clumped.fastq in2=${accno}_2.clumped.fastq out=${accno}_1.clumped.trimmed.fastq out2=${accno}_2.clumped.trimmed.fastq ziplevel=5 pigz ordered qtrim=rl trimq=15 minlen=15 ecco=t maxns=5 trimpolya=10 1>>pipeline_$accno.log 2>&1


echo "" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "       + BBTools: BBNorm (Coverage Normalization) Started                       +" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "" >> pipeline_$accno.log

$toolbox_path/bbmap/bbnorm.sh in=${accno}_1.clumped.trimmed.fastq in2=${accno}_2.clumped.trimmed.fastq out=${accno}_1.clumped.trimmed.normalized.fastq out2=${accno}_2.clumped.trimmed.normalized.fastq target=10 min=2 histcol=2 khist=${accno}.clumped.trimmed.normalized.khist peaks=${accno}.clumped.trimmed.normalized.peaks 1>>pipeline_$accno.log 2>&1


echo "" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "       + Lighter: Lighter (Error Correction) Started                            +" >> pipeline_$accno.log
echo "       --------------------------------------------------------------------------" >> pipeline_$accno.log
echo "" >> pipeline_$accno.log

$toolbox_path/lighter/lighter -r ${accno}_1.clumped.trimmed.normalized.fastq -r ${accno}_2.clumped.trimmed.normalized.fastq -K 27 100000000 -t 4 -zlib 5


echo "" >> pipeline_$accno.log
echo "########################################################################################" >> pipeline_$accno.log
echo "# Preprocessing Pipeline Finished                                                      #" >> pipeline_$accno.log
echo "########################################################################################" >> pipeline_$accno.log