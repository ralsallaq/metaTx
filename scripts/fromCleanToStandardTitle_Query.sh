#!/bin/bash
#requires python
module load python/2.7.13
nFs=`ls -p alignedQuery/sample_*.db.qc.clean.align.fasta | wc -l`
#$1=refalignment file with clean seq titles acceptable vy RAxML, e.g. dbA.taxa.rax.clean.align.fasta
#$2=the original reference file with ncbi titles , e.g. dbA.taxa.fasta
dir=`echo $PWD`
function func1 {
    local i=$1
    file=`ls -p alignedQuery/sample_*.db.qc.clean.align.fasta | head -$i |tail -1`
    prefix=`echo $file |cut -f 2 -d "_"|cut -f 1 -d "."`
    echo "processing $prefix"
    if [ ! -f raxml_heads ]; then
        grep "^>" $1 |cut -f2 -d ">" > raxml_heads
    fi
    if [ ! -f ncbi_heads ]; then
        grep "^>" $2 |cut -f2 -d ">" >ncbi_heads
    fi
    save=`echo $file|sed "s/\.fasta/\.fasta_save/g"`
    cp -f $file $save
    out=`echo $file|sed "s/rax\.//g;s/\.fasta/\.ncbi\.fasta/g"`
    sed '/^>/s/_description_.*organism/_organism/g' $file > temp_$prefix
    python $dir/scripts/fromCleanToStandardTitle.py -i1 raxml_heads -i2 ncbi_heads -iF temp_$prefix -oF $out
}
    
export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -n 1 -R "span[hosts=1] rusage[mem=2000]" -P microbiome -J "align[1-$nFs]" -oo "logging/cleanAlignQuery_bysample.log" "func1 \$LSB_JOBINDEX" 

#for (( i=1;i<=$nFs;i++ )); do
#    func1 $i
#done
