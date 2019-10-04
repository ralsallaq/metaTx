#!/bin/bash
#requires python 
module load python/2.7.13
if [ ! -d "analysis" ]; then
    mkdir analysis
fi
files=`echo $@ | tr " " "\n"|grep '.*.jplace'`
outFs=`echo $@ | tr " " "\n"|grep '.*.csv'`
#number of sample files
nFs=`echo $files|wc -w`
echo $nFs
echo $files

##go through the sample files
function func1 {
    local j=$1
    local fls=$2
    file=`echo $fls | cut -f $j -d " "`
    outPrefix=`basename $file .jplace|sed -e 's/sample_//g'`

    echo "save a csv file for pplacer statistics from jplace file for $outPrefix"
    echo  "the command is python ./scripts/jplaceStats.py $file analysis/'$outPrefix'_pplaceStats.csv"

    python ./scripts/jplaceStats.py $file --csv analysis/"$outPrefix"_pplaceStats.csv

    echo $outPrefix

}

export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=4000]"  -P microbiome -J "jplaceStats[1-$nFs]" -oo "logging/jplaceStats_%I.log" "func1 \$LSB_JOBINDEX '$files'"
