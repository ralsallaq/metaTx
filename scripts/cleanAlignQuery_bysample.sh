#!/bin/bash
#requires seqmagick and infernal 
module load python/2.7.13
module load infernal/1.1.2
module load fasttree/2.1.10
if [ ! -d "alignedQuery" ]; then
    mkdir alignedQuery
fi
if [ ! -f 16S_bacteria.cm ]; then
    echo "Missing covariance model file 16S_bacteria.cm"
fi
raxrefAlign=$1
#shift to the next argument
shift
#all the files to be processed
files=$@
#number of sample files
nFs=`echo $files |wc -w`
echo $nFs

#go through the sample files
function func1 {
    local i=$1
    local files=$2
    local raxrefAlign=$3
    file=`echo $files|cut -f $i -d " "`
    outPrefix=`basename $file|cut -f 1 -d "."`

    echo "------processing file : $file-----------------"
    echo $outPrefix
    echo "----quality filtering ------"
    echo "----cleaning punctuations from $file------------"
    ## raxml Illegal characters in taxon-names are: tabulators, carriage returns, spaces, ":", ",", ")", "(", ";", "]", "[", "'"

    tr "[:\ \' \%\* + (),;\\ \"{}=]" "-" < $file > temp_$outPrefix 
    tr -d "_" < temp_$outPrefix > temp2_$outPrefix
    tr "\." "-" < temp2_$outPrefix > temp3_$outPrefix
    sed '/^>/s/-/_/g;/^>/s/___/_/g;/^>/s/____/_/g;/^>/s/_$//g;/^>/s/__/_/g' temp3_$outPrefix > alignedQuery/$outPrefix.assembled.qc.clean.fasta
    wait
    rm -f temp_$outPrefix temp2_$outPrefix temp3_$outPrefix
    echo "---for each sample concatentate query reads with the sample refrence alinged sequences-----------" 
    cat $raxrefAlign  alignedQuery/$outPrefix.assembled.qc.clean.fasta > alignedQuery/sample_$outPrefix.db.qc.clean.fasta 
    ##remove gaps in place using sqmagick
    seqmagick mogrify --ungap alignedQuery/sample_$outPrefix.db.qc.clean.fasta
    wait
    echo "------aligning query+ref using infernal------------"
    cmalign --mxsize 9000  --cpu 4 --dnaout -o alignedQuery/sample_$outPrefix.db.qc.clean.align.sto --outformat Pfam 16S_bacteria.cm alignedQuery/sample_$outPrefix.db.qc.clean.fasta >alignedQuery/cmalingout_$outPrefix.log

    wait
    echo "----------------converting to fasta -----------"
    if [ -f alignedQuery/sample_$outPrefix.db.qc.clean.align.sto ]; then 
     #convert to fasta
     seqmagick convert alignedQuery/sample_$outPrefix.db.qc.clean.align.sto alignedQuery/sample_$outPrefix.db.qc.clean.align.fasta 
    fi

    echo $outPrefix

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard"  -n 4 -R "span[hosts=1] rusage[mem=12000]" -P microbiome -J "align[1-$nFs]" -oo "logging/alignedQuery_%I.log" "func1 \$LSB_JOBINDEX '$files' '$raxrefAlign' " 
