#!/bin/bash
#requires: raxml
module load raxml/AVX/8.0.22

alignedcleanAforRaXML=$1

function func1 {
    local  alignedcleanAforRaXML=$1 
    if [ ! -f $alignedcleanAforRaXML ]; then
        echo "are you sure you have a clean aligned ref fasta file?edy for raxml?"
    fi
    echo "------Building a tree without confidence scores just to get stats for Taxtastic------------"
    echo "----- look into the RAxML*info*tre file and how many distinct patterns are in their: choose one core/thread per 500 DNA patterns---"
    echo "----- see page 5 of this manual https://sco.h-its.org/exelixis/resource/download/NewManual.pdf---"

    if [ -f $alignedcleanAforRaXML ]; then

        raxmlHPC-PTHREADS-AVX -T 8 -m GTRGAMMA -s  $alignedcleanAforRaXML  -n ref.tre -f d -p 12345
        wait
        if [ -f "RAxML_bestTree.ref.tre" ]; then
            echo "-----Root the best tree obtained from previous step using -f I method---------------------------"
            #raxmlHPC-PTHREADS-AVX -T 8 -m GTRGAMMA -f I -t RAxML_bestTree.ref.tre -n root.ref.tre 
            raxmlHPC-PTHREADS-AVX -T 8 -m GTRGAMMA -f I -t RAxML_bestTree.ref.tre -n root.ref.tre 
        fi
        wait
        if [ -f "RAxML_rootedTree.root.ref.tre" ]; then
            echo "------Generate confidence scores for the rooted tree ---------------"
            raxmlHPC-PTHREADS-AVX -T 8 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root.ref.tre -n conf.root.ref.tre -s $alignedcleanAforRaXML
        fi
    fi
    echo "----now you have all what you need to build a reference package and then start processing query reads ---------"

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "rhel7_large_mem" -n 8 -R "span[hosts=1] rusage[mem=20000]" -P microbiome -J "raxml" -oo "logging/buildRefTree.log"  "func1 $alignedcleanAforRaXML"
