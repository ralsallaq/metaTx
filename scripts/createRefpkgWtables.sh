#!/bin/bash
#requires: taxtastic 
module load python/2.7.13

#$1 seqInfo file dbA.taxa.SeqInfo.csv
#$2 taxonomyInfo file dbA.taxta.TaxaInfo.csv
#$3 tree file RAxML_fastTreeSH_Support.conf.root.ref.nwk
#$4 tree stats file RAxML_info.ref.tre
#$5 alignment ref fasta with primary taxids dbA.taxa.ncbi.clean.align.primaryTaxids.fasta

function func1 {
   
    local seqInfo=$1
    local taxInfo=$2
    local tree=$3
    local stats=$4
    local refalgn=$5


    echo "----use the reference alignment, the rooted tree with confidence scores for reference alignment, and the raxml stats file without the confidence scores to create the reference package that will be used for all samples ------------------" 
    echo "to use taxonomy tables for the reference alignments instead of their original titles you should use the following command: look here to see how seqinfo and taxatable look like https://www.dropbox.com/sh/lqj7i2ra46s6fsu/AADzY8j05NhonARLNrFZ7A-pa?dl=0"
    echo "to get the taxonomy table see https://groups.google.com/forum/#!searchin/pplacer-users/taxonomy$20file%7Csort:date/pplacer-users/5TUDgLSuozs/aCC0u2XCAQAJ"

    echo "check all good"
    seqInls=`wc -l $seqInfo | cut -f1 -d " "`
    refalgnnls=`grep '^>' $refalgn|wc -l`
    taxInfonls=`wc -l $taxInfo | cut -f1 -d " "`
    if [[ $((seqInls-1)) -ne $((refalgnnls)) ]]; then
        echo "seems the sequence information file is not complete ...check"
        exit 1
    fi
    if [[ $((taxInfonls)) -eq 0 ]]; then
        echo "the taxinfo file is empty ..check"
        exit 1
    fi



    taxit create -l 16S_rRNA -P allsamples.refpkg --taxonomy $taxInfo --seq-info $seqInfo --aln-fasta $refalgn --tree-stats $stats --tree-file $tree

    echo "---done with the reference package, now carry one with processing samples ------"
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=2500]" -P microbiome -J "taxit" -oo "logging/createRefpkgWtables.log" "func1 $1 $2 $3 $4 $5"
