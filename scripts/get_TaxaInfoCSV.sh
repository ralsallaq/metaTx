#!/bin/bash
#requires taxtastic
#activate virtual env for taxtastic
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
source $activate taxtastic
#echo "Do not forget to activate taxtastic environment on conda using: conda activate taxtastic"
#$1 is the seqInfo file: e.g dbA.taxa.SeqInfo.csv
#$2 taxonomy database
#$3 is the output TaxInfo file

function func1 {
    local    seqInfo=$1
    local    taxDB=$2
    local    taxInfo=$3
    if [ ! -f $taxInfo ]; then
        if [ -f $taxDB ]; then
            awk -F "\"*,\"*" '{print $2}' $seqInfo > temp_ids
            sed 's/tax_id//g' temp_ids > temp_ids_
            sed '/^$/d' temp_ids_ > taxids
            rm -f temp_ids_ temp_ids
            taxit taxtable $taxDB -f taxids -o $taxInfo
        
        else
            echo "taxonomy.db not found download and save via taxtastic using:"
            echo "taxit new_database taxonomy.db "
        fi
    else
        echo "file $taxInfo already exists"
    fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=6000]" -P microbiome -J "TaxInfo" -oo "logging/get_TaxaInfoCSV.log" "func1 $1 $2 $3"
