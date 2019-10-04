#!/bin/bash
#requires R for dada2 and dada2
#load R by load python/2.7.15-rhel7
module load python/2.7.15-rhel7
#$1 is the absolute path to the R script for applying dada2 analysis
#$2 is the absolute path to the sampleInfo CSV file with Path, R1 R2  dirName and Sample_Name fields   
#$3 is the absolute path to the fasta files encompassing SVs for each sample
function func1 {
local dada2RScript=$1
local sampleInfoCSV=$2
local SVdir=$3

if [ ! -f $sampleInfoCSV ]; then
    echo "The sample information CSV file $1 you specified  does not exist please check its path!"
    exit 1
fi
if [ ! -d $SVdir ]; then
    Rscript $dada2RScript $sampleInfoCSV $SVdir
else
    echo "the SVs directory already exist, doing nothing"
fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1 

bsub -q "standard" -n 8 -R "span[hosts=1] rusage[mem=8000]" -P microbiome -J "dada2" -oo "logging/get_dada2_svs.log" "func1 $1 $2 $3"
