#!/bin/bash
#requires pplacer/guppy/rppr
module load pplacer
#This bash script executes the phylogenetic placement using pplacer  
if [ ! -d "phyloPlace" ]; then
    echo "seems you have not done the pplacer step"
fi

if [ ! -d "analysis" ]; then
    mkdir analysis
fi

Ncpu=4
refpkgdir=`echo $@| tr " " "\n"|grep '.*.refpkg'`
algnFs=`echo $@ | tr " " "\n"|grep '.*.fasta'`
files=`echo $@ | tr " " "\n"|grep '.*.jplace'`
#number of sample files
nFs=`echo $files|wc -w`

echo $nFs

##go through the sample files
function func1 {

    i=$1
    files=$2
    file=`echo $files| cut -d " " -f $i`
    outPrefix=`basename $file .jplace|sed -e 's/sample_//g'`
    refpkgdir=$3
    algnFs=$4
    #algnF=`echo $algnFs| tr " " "\n"|grep "$outPrefix"`
    algnF=`echo alignedQuery/sample_$outPrefix.db.qc.clean.align.ncbi.fasta`
    Ncpu=$5


    echo $refpkgdir $algnFs $files 
    echo "processing these two files---"
    echo $algnF $file
    echo "-----------"
    echo "------processing file : $file-----------------"
    echo $outPrefix
#    echo "Generate an easily parsed csv file of placements, with with only a single placement reported for each query read"
#    if [ ! -f analysis/$outPrefix.csv ]; then
#        guppy to_csv -o analysis/$outPrefix.csv --pp --point-mass $file 
#    fi
#    echo "Generate a phyloxml tree with edges fattened according to the number of placements"
#    if [ ! -f analysis/$outPrefix.phyloxml ]; then
#        guppy fat -o analysis/$outPrefix.phyloxml --node-numbers --point-mass --pp $file 
#    fi
    echo "Generate a csv for the edpl (expected distance b/w placement locations metric which summarizes the uncertainty in placement using pplacer; it is extremely useful when there are a number of closely-related sequences in the reference alignment (subspecies or strains); it measures the average resident time on an edge (the higher the metric the higher the time)"
    if [ ! -f analysis/$outPrefix.edpl.csv ];then
        guppy edpl  -o analysis/$outPrefix.edpl.csv --csv --pp --first-only $file 
    fi
    echo "Generate a sqlite database of placements"
    #using hybrid2 classifier which combines Neighborhood-Based-Clustering NBC and tree-based algorithm 
    echo "first prepare the sqlite3 database to save placement summaries populate taxa and ranks from the refpkg"
    if [ ! -f analysis/$outPrefix.db ]; then
        rppr prep_db -c $refpkgdir --sqlite analysis/$outPrefix.db
        wait
        echo "then save other tables into the database regarding placements"
        guppy classify --pp -j $Ncpu --classifier hybrid2 --nbc-sequences $algnF  -c $refpkgdir --sqlite analysis/$outPrefix.db $file >analysis/out_$outPrefix.classify
    fi

    echo "calculate the ADCL metric which is the average distance to closest leaf"
    if [ ! -f analysis/$outPrefix.adcl.csv ]; then
        guppy adcl -o analysis/$outPrefix.adcl.csv $file --pp --no-collapse >analysis/out_$outPrefix.adcl.out
    fi

    echo $outPrefix

}

export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 4 -R "span[hosts=1] rusage[mem=10000]"  -P microbiome -J "phyloPlace_classify[1-$nFs]" -o "logging/phyloPlace_classify.%I" "func1 \$LSB_JOBINDEX '$files' '$refpkgdir' '$algnFs' '$Ncpu'"

