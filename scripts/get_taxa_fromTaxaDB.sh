#!/bin/bash
#requires python
module load python/3.5.2
taxaout=$1
taxTable=$2
taxDB=$3

bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "getTaxa" -oo "logging/get_taxa_fromTaxaDB.log" "python $PWD/scripts/get_taxa_fromTaxaDB.py -db $taxDB -l $taxTable -o $taxaout"
