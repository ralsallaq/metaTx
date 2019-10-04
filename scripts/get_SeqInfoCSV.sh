#!/bin/bash
#requires python
module load python/2.7.13
#$1 is the ref fasta file with ncbi titles e.g. dbA.taxa.ncbi.clean.align.fasta
#$2 email for NCBI server
#$3 is the name of the seqInfo CSV file
#$out is the ref fasta file with ncbi titles and taxids changed to primary taxids (if any)
dir=`echo $PWD`
baseN=`basename $1 .fasta`
save=`echo $baseN".fasta_save"`
cp -f $1 $save
out=`echo $baseN".primaryTaxids.fasta"`
if [ ! -f $3 ]; then
    bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "SeqInfo" -oo "logging/get_SeqInfoCSV.log" "python $dir/scripts/get_seqInfo_ncbi.py -i $1 -e $2 -oC $3 -oF $out"
else
    echo "file $3 already exists"
fi
