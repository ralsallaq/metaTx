#!/bin/bash
#bash script
files=`ls -p analysis/*.db`

function func1 {
local files=$1
nfiles=`echo $files|wc -w` 
arbitfile=`ls -p analysis/*.db|head -3|tail -1`
dbName=`sqlite3 -header -column $arbitfile "select * from multiclass"|awk '{print $2}'|tail -1`
#echo $dbName
res=`grep "^>$dbName" alignedQuery/*.assembled.qc.clean.fasta`
if [ "$res" ]; then
    x=0
else
    echo "check if $arbitfile has the same seqname as in alignedQuery/*.assembled.qc.clean.fasta; because they seem to be different!" 
    x=1
    exit 1
fi
echo "number of files $nfiles"
for ((i=1;i<=$nfiles;i++)); do
    f=`echo $files|cut -f $i -d " "`
    nlines=`sqlite3 -header -column $f "select * from multiclass"|wc -l`
    if [[ $nlines -lt 1 ]]; then
        echo "dbase file $f is empty...check"
        x=1
        exit 1 
        break
    fi
done
if [[ $x = 0 ]]; then
echo "check point passed successfully"
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "checkdbnames" -oo "logging/chkdbnames.log" "func1 '$files'"
