#!/bin/bash
#requirs: seqmagick and python
module load python/2.7.13 
#dbvs=$1
#dbfull=$2 #either a fasta file or a seqInfo CSV file
#out=$3
function func1 {
local dbvs=$1
local dbfull=$2
local out=$3
if [ ! -f $out ]; then
    nFields=`cat "$dbfull" | awk 'BEGIN{FS=","}END{print NF}'`
    nFastaL=`grep '^>' "$dbfull"|wc -l`
    if [[ "$nFields" -gt 1 ]]; then
        echo "Using CSV seqInfo file to assign taxonomy"
        cmd=./scripts/assignTaxonomyToVsearchFasta_useCSVF.py
    elif [[ "$nFastaL" -gt 0 ]]; then
        echo "Using NCBI formatted fasta file to assign taxonomy"
        cmd=./scripts/assignTaxonomyToVsearchFasta_useFastaF.py
    else
        echo "the file that is passed to assign taxa is neither a csv file nor a fasta file"
        exit 1
    fi

    python $cmd -v $dbvs -f $dbfull -o temp.fasta 
    echo "remove double quotes around organism name"
    sed -e 's/organism=""/organism="/g; s/"" taxonomy/" taxonomy/g' temp.fasta>$out
    wait
    echo "dereplicate in place any duplicated sequences"
    seqmagick mogrify --deduplicate-sequences $out
    echo "done"
else
    echo "file $out already exists"
fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=6000]" -P microbiome -J "assignTaxa" -oo "logging/assignTaxaToVsearch.log" "func1 $1 $2 $3"
