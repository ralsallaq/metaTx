#!/bin/bash
#requires python 
module load python/3.5.2
dir=`echo $PWD`
saveFile=$1
svSeqTab=$2
shift
shift
files=$@
#number of sample files
nFs=`echo $files|wc -w`
echo $nFs
if [ ! -f $saveFile ]; then
    bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=6000]"  -P microbiome -J "oneTaxa" -oo "logging/tallySamplesByTaxa_keepBest.log" "python $dir/scripts/tallySamplesByTaxa_keepBest.py '$saveFile' '$svSeqTab' '$files'"
else
    echo "file $saveFile already exists"
fi
