#!/bin/bash
#bash script
function func1 {
alignedcleanA=$1
baseN=`basename $alignedcleanA .fasta`
alignedcleanAforRaXML=`echo $baseN".rax.fasta"`
echo "------------------Clean alignment file for raxml--------"
sed '/^>/s/_description_.*$//g;/^>/s/_organism_.*//g' $alignedcleanA >temp_fasta
wait
mv -f temp_fasta  $alignedcleanAforRaXML 
}

export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "prepForRax" -oo "logging/prepForRax.log" "func1 $1"
