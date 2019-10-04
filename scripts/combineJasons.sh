#!/bin/bash
#requires pplacer/guppy
module load pplacer/1.1.alpha19
if [ ! -d "phyloPlace" ]; then
    echo "seems you have not done the pplacer step, ...stopping"
    exit 1
fi
if [ ! -d "analysis" ]; then
    mkdir analysis
fi
combinedJason=$1
prefix=`basename combinedJason .jplace`
shift  
files=$@ #`echo $@ | tr " " "\n"|grep '*.jplace'`
#number of sample files
nFs=`echo $files|wc -w`
echo "number of jplace files to combine = $nFs"

function func1 {

    combFile=$1
    files=$2
    prefix=$3

    echo $combFile $files 
    echo "-----------"

    echo "merge jplace files into one and produce a csv file indicating the source of each read in the combined placefile"

    if [ ! -f $combFile ]; then
        guppy merge -o $combFile --split-csv  $files
    else
        echo "file $combFile already exists"
    fi

}

export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 4 -R "span[hosts=1] rusage[mem=6000]"  -P microbiome -J "combineJplace" -oo "logging/combineJason.log" "func1  '$combinedJason' '$files' '$prefix'"
