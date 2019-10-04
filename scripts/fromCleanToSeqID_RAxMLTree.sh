#!/bin/bash
#requires python
module load python/2.7.13
#$1=raxml newick tree with support info and tip-labels acceptable by RAxML, e.g. RAxML_fastTreeSH_Support.conf.root.ref.tre
#$2=raxml_heads
#$3=ncbi_heads
#$out=newick tree with support info and tip-labels changed to seqID as in the SeqInfo csv file, e.g. RAxML_fastTreeSH_Support.conf.root.ref.nwk
if [ -f RAxML_fastTreeSH_Support.conf.root.ref.nwk ]; then
    echo "RAxML_fastTreeSH_Support.conf.root.ref.nwk exists. Stop."
    exit 1
fi
dir=`echo $PWD`
save=`echo $1|sed "s/ref\.tre/ref\.tre_save/g"`
cp -f $1 $save
out=`echo $1|sed "s/ref\.tre/ref\.nwk/g"`
bsub -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "fixTipLabels" -oo "logging/fromCleanToSeqID_RAxMLTree.log" "python $dir/scripts/fromCleanToSeqID_RAxMLTree.py -i1 $2 -i2 $3 -iT $1 -oT $out"
