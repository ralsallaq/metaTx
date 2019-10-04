#!/bin/bash
#requires vsearch 
module load vsearch/2.8.0
masked_reffasta=$1 #database
#shift to the next argument
shift
#all the files to be processed
files=$@
#number of sample files
nFs=`echo $files |wc -w`
echo $nFs

if [ ! -d "vsearchOut" ]; then
    mkdir vsearchOut
fi

function func1 {
    local i=$1
    local files=$2
    local masked_reffasta=$3
    file=`echo $files|cut -f $i -d " "`
    #this part assumes that dada2 fasta are named as sampleName.fasta
    outPrefix=`basename $file|cut -f 1 -d "."`
    echo "------processing sample $outPrefix in file: $file-----------------"
    #Search for matches with a threshold here I am using a low threshold of 80%, maxaccepts and maxrejects is set to 0 to consider all database sequences
    vsearch --threads 4 --maxaccepts 20 --usearch_global $file --db $masked_reffasta --id 0.8 --dbmatched vsearchOut/matched_$outPrefix.fasta --alnout vsearchOut/matched_$outPrefix.txt 
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard" -n 4 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "vsearch[1-$nFs]" -oo "logging/vsearch_%I.log" "func1 \$LSB_JOBINDEX '$files' '$masked_reffasta'"
