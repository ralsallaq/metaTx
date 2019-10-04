#!/bin/bash
#requires seqmagick and infernal  
module load python/2.7.13
module load infernal/1.1.2
module load fasttree/2.1.10
refdbFasta=$1
if [ ! -f 16S_bacteria.cm ]; then
    echo "Missing covariance model file 16S_bacteria.cm ... stopping"
    exit 1
fi
if [ ! -f $refdbFasta ]; then
    echo " missing the file $refdbFasta ... stopping"
    exit 1
fi
if [ -f "dbvs.clean.align.sto" ]; then
    echo "seems you already have an alignment reference database...if this is a mistake delete the dbA.clean.align.sto file and re-run"
    echo "..........stopping"
    exit 1
else
    function func1 {
        local refdbFasta=$1
        baseN=`basename $refdbFasta .fasta`
        cleanF=`echo $baseN".clean.fasta"` 
        cleanalgnS=`echo $baseN".clean.align.sto"`
        cleanalgnF=`echo $baseN".clean.align.fasta"`
        echo "------processing file : $refdbFasta-----------------"
        echo "----cleaning punctuations------------"
        tr "[:\ \' \%\* + (),;\\ \"{}=]" "-" < $refdbFasta > temp_fasta 
        tr -d "_" < temp_fasta > temp2_fasta
        tr "\." "-" < temp2_fasta > temp3_fasta
        sed '/^>/s/-/_/g;/^>/s/___/_/g;/^>/s/____/_/g;/^>/s/_$//g;/^>/s/__/_/g;/^>/s/_description_organism//g' temp3_fasta > $cleanF 
        wait
        rm -f temp_fasta temp2_fasta  temp3_fasta
        echo "------aligning using infernal------------"
        cmalign  --mxsize 9000 --cpu 8 --dnaout -o $cleanalgnS --outformat Pfam 16S_bacteria.cm $cleanF > cmalignRef.out
        wait
        echo "----------------converting to fasta -----------"
        if [ -f $cleanalgnS ]; then 
         #convert to fasta
         seqmagick convert $cleanalgnS $cleanalgnF 
        fi
        wait
        echo "   ------deduplicate inplace using seqmagick------------"
        if [ -f $cleanalgnF ]; then
            seqmagick mogrify --deduplicate-sequences dbA.taxa.clean.align.fasta
        fi
        echo "-----done with: now you have a deduplicated reference alignment fasta $cleanalgnF-----------------"

    }
    export LSB_JOB_REPORT_MAIL="N"
    export -f func1

    bsub -q "rhel7_large_mem" -n 8 -R "span[hosts=1] rusage[mem=18000]" -P microbiome -J "alignment" -oo "logging/cleanAlignRef.log" "func1 '$refdbFasta'"
fi
