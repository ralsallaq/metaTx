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
echo "input files are  $saveFile $svSeqTab $files" 
#number of sample files
echo "$nFs"
if [ ! -f $saveFile ]; then 
    bsub -n 1 -R "span[hosts=1] rusage[mem=2000]"  -P microbiome -J "multipleTaxa" -oo "logging/tallySamplesByTaxa.log" "python $dir/scripts/tallySamplesByTaxa.py '$saveFile' '$svSeqTab' '$files'"
else
    echo "file $saveFile already exists" 
fi
