#!/bin/bash
#requires python
module load python/2.7.13
#$1=refalignment file with clean seq titles acceptable vy RAxML, e.g. dbA.taxa.rax.clean.align.fasta
#$2=the original reference file with ncbi titles , e.g. dbA.taxa.fasta
dir=`echo $PWD`

if [ ! -f raxml_heads ]; then
    grep "^>" $1 |cut -f2 -d ">" > raxml_heads
elif [ ! `wc -l raxml_heads|cut -f 1 -d " "` -gt 0 ]; then
    grep "^>" $1 |cut -f2 -d ">" > raxml_heads
fi
if [ ! -f ncbi_heads ]; then
    grep "^>" $2 |cut -f2 -d ">" >ncbi_heads
elif [ ! `wc -l ncbi_heads|cut -f 1 -d " "` -gt 0 ]; then
    grep "^>" $2 |cut -f2 -d ">" >ncbi_heads
fi


nlrax=`wc -l raxml_heads|cut -f 1 -d " "`
nlncbi=`wc -l ncbi_heads|cut -f 1 -d " "`
if [ ! $nlrax -gt 0 ]; then
    grep "^>" $1 |cut -f2 -d ">" > raxml_heads
fi
if [ ! $nlncbi -gt 0 ]; then
    grep "^>" $2 |cut -f2 -d ">" >ncbi_heads
fi
save=`echo $1|sed "s/\.fasta/\.fasta_save/g"`
cp -f $1 $save
out=`echo $1|sed "s/rax\.//g;s/\.fasta/\.ncbi\.fasta/g"`
if [ ! -f "$out" ]; then
    bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "fixTitles" -oo "logging/fromCleanToStandardTitle.log" "python $dir/scripts/fromCleanToStandardTitle.py -i1 raxml_heads -i2 ncbi_heads -iF $1 -oF $out"
else
    echo "file $out already exists"
fi
