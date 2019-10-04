#!/bin/bash
#requires: python
module load python/2.7.13
files=`echo $@ | tr " " "\n"|grep '.*.db'`
qureyFs=`echo $@ | tr " " "\n"|grep '.*.fasta'`
#number of sample files
nFs=`echo $files|wc -w`
echo $nFs
echo $files

##go through the sample files
function func1 {
    local j=$1
    local fls=$2
    local qfls=$3
    file=`echo $fls | cut -f $j -d " "`
    echo $file $queryF
    outPrefix=`basename $file .db`
    queryF=`echo alignedQuery/$outPrefix.assembled.qc.clean.fasta`

    echo "save a csv file for the best rank and its statistics from db for $outPrefix"
    echo  "the command is python ./scripts/bestRankForReads.py -db $file -oc analysis/'$outPrefix'_bestRankStats.csv"

    python ./scripts/bestRankForReads.py -db $file -q $queryF -oc analysis/"$outPrefix"_bestRankStats.csv

    echo $outPrefix

}

export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=4000]"  -P microbiome -J "bestRankStats[1-$nFs]" -oo "logging/bestRankStats_%I.log" "func1 \$LSB_JOBINDEX '$files' '$qureyFs'"
#for ((i=1 ; i<=$nFs ; i++)) ; do
#    func1 "$i" "$files"
#done
