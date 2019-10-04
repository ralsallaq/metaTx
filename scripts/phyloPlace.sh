#!/bin/bash
#requires python and pplacer
#This bash script invokes pplacer to place aligned query reads (which is aligned with the reference sequences) onto a reference package 
module load python/2.7.13
module load pplacer
if [ ! -d "phyloPlace" ]; then
    mkdir phyloPlace
fi
#first argument
refpkgdir=$1
shift
#everything else
files=$@
#number of sample files
nFs=`echo $files|wc -w`
echo $nFs
##go through the sample files

function func1 {
    local i=$1
    local files=$2
    local refpkgdir=$3
    dir=`echo $PWD`
    file=`echo $files|cut -f $i -d " "`
    outPrefix=`basename $file|cut -f 1 -d "."|sed -e 's/sample_//g'`
    #outPrefix=`echo $file |cut -f 2 -d "/"|cut -f 1 -d "."|cut -f 2 -d "_"`
    if [ ! -f phyloPlace/sample_"$outPrefix".jplace ]; then
        echo "------processing file : $file-----------------"
        echo $outPrefix
    
        echo "-----first convert titles in the alignment query file to standard ncbi titles-------"
        save=`echo $file|sed "s/\.fasta/\.fasta_save/g"`
        cp -f $file $save
        out1=`echo $file|sed "s/rax\.//g;s/\.fasta/\.ncbi\.fasta/g"`
        out2=`echo $out1|sed "s/\.fasta/\.primaryTaxid\.fasta/g"`
        #sed '/^>/s/_description_.*organism/_organism/g' $file > temp_$outPrefix
        python $dir/scripts/fromCleanToStandardTitle.py -i1 raxml_heads -i2 ncbi_heads -iF $file -oF $out1
        wait
    
    
        echo "-----second change any bad taxids into good taxids -----"
        if [ -f nomatchedTaxIDs.csv ]; then
            echo "there is some bad taxids with information in nomatchedTaxIDs.csv; titles of ref sequences  in query alignments will be changed to primary taxids"
            nmismatches=`wc -l nomatchedTaxIDs.csv|cut -f 1 -d " "`
            #replace bad taxids with good taxids sequentially 
            cp $out1 temp_$outPrefix
            wait
            for ((m=1;m<=$((nmismatches-1));m++)); do
                cp temp_$outPrefix temp1_$outPrefix
                badTaxId=`awk -F "\"*,\"*" '{print $1}' nomatchedTaxIDs.csv|sed 's/tax_id//g;s/^"//g;/^$/d'|head -$m|tail -1`
                goodTaxId=`awk -F "\"*,\"*" '{print $3}' nomatchedTaxIDs.csv|sed 's/primary_tax_id//g;s/^"//g;/^$/d'|head -$m|tail -1`
                wait
                sed 's/ncbi_tax_id=\"'"$badTaxId"'\"/ncbi_tax_id=\"'"$goodTaxId"'\"/g' temp1_$outPrefix > temp_$outPrefix 
            done
            cp temp_$outPrefix $out2
            wait
            pplacer --timing -j 4  --prior-lower 0.01 --inform-prior --map-identity -o phyloPlace/sample_$outPrefix.jplace -p --keep-at-most 20 -c $refpkgdir $out2 > phyloPlace/pplacerout_$outPrefix
        else
            echo "all taxids are primary taxids and no change in the alignment query files is necessary"
            pplacer --timing -j 4  --prior-lower 0.01 --inform-prior --map-identity -o phyloPlace/sample_$outPrefix.jplace -p --keep-at-most 20 -c $refpkgdir $out1 > phyloPlace/pplacerout_$outPrefix
        fi
    
        echo $outPrefix
    else
        echo "file phyloPlace/sample_$outPrefix.jplace already exists" 
    fi


}

export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 4 -R "span[hosts=1] rusage[mem=12000]"  -P microbiome -J "phyloPlace[1-$nFs]" -oo "logging/phyloPlace\%I.log" "func1 \$LSB_JOBINDEX '$files' '$refpkgdir'"

