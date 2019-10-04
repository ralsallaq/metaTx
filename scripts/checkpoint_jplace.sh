#!/bin/bash
#bash script 
files=`ls -p phyloPlace/sample_*.jplace`

function func1 {
local files=$1
x=0
for f in $files; do 
sampleID=`basename $f '.jplace'|sed -e 's/sample_//g'`
if [ ! -f alignedQuery/sample_$sampleID.db.qc.clean.align.fasta ]; then
    x=1
    echo "missing file alignedQuery/sample_$sampleID.db.qc.clean.align.fasta"
fi
checkFile=`grep "tree" $f|wc -l`
if [[ $checkFile -lt 1 ]]; then
    x=1
    echo "check file $f, it seems it is erroneous" 
fi
done

if [[ $x = 0 ]]; then
    echo "checkpoint successfully passed"
else
    echo "checkpoint failed, it seems the jplace files are named differently than the alignment files...check"
    echo "...or some of the place files are erroneous"
    exit 1
fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "checkjplace" -oo "logging/chkjplace.log" "func1 '$files'"
