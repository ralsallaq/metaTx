#!/bin/bash
#requires pplacer
module load pplacer
refpkgdir=`echo $@| tr " " "\n"|grep '.*.refpkg'`
files=`echo $@ | tr " " "\n"|grep '.*.jplace'`
#number of sample files
nFs=`echo $files|wc -w`

echo $nFs

##go through the sample files
function func1 {

    i=$1
    files=$2
    file=`echo $files| cut -d " " -f $i`
    outPrefix=`basename $file .jplace|cut -f 2 -d "_"`
    refpkgdir=$3

    echo $refpkgdir $files $file 
    echo "------processing file : $file-----------------"
    echo $outPrefix
    echo "find a good collection of sequences to cut from a placefile's ref tree"
    rppr min_adcl --point-mass --pp -c $refpkgdir -o  analysis/$outPrefix.min_adcl.csv -t analysis/$outPrefix.trimmedTreeByMinAdcl.tre --max-adcl 0.0  $file >analysis/$outPrefix.min_adcl.out
    echo $outPrefix

}

export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -n 1 -R "span[hosts=1] rusage[mem=10000]"  -P microbiome -J "minADCL[1-$nFs]" -oo "logging/get_minADCL.log" "func1 \$LSB_JOBINDEX '$files' '$refpkgdir'"
