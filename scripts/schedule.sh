#!/bin/bash
#bash script to schedule the pipeline on the cluster using gnu-make and lsf (IBM)

usage()
{
    echo "usage:./scripts/schedule.sh -r bigDB -f querySampleInfo -d taxDB -e myemail"

} 


while getopts ":r: :f: :d: :e:" opt; do
    case $opt in
        r)
            echo "The large reference fasta file is set to $OPTARG"
            bigDB=$OPTARG
            ;;
        f)
            echo "The sample info csv file is set to $OPTARG"
            querySampleInfoF=$OPTARG
            ;;
        d)
            echo "The sqlite3 taxonomy database is set to $OPTARG"
            taxDB=$OPTARG
            ;;
        e)
            echo "The email for NCBI server inquiries is set to $OPTARG"
            myemail=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG reguires an argument"
            exit 1
            ;;
    esac
done

#test if options are passed
if [ -z "$bigDB" ] || [ -z "$querySampleInfoF" ] || [ -z "$taxDB" ]; then
    usage
    exit 1
fi


function checkJoboutput {
local output=$1
echo "This job output is $output"
check=`echo $output | grep 'make\|Job' |awk '{print $NF}'`
echo check is $check
{ # previous command failed try this next: check if make returns up to date
    #check=`echo $output | awk '{print $NF}'`
    #echo "check is $check"
    if [[ $check == "date." ]]; then
        echo "The command returned up to date status"
        continue
    elif [[ -z $check ]]; then 
        echo "The command failed"
        exit 1
    else
        false
    fi
} || \
{ #try to see if it is a scheduled job array on the cluster
    if [[ -z $check ]]; then
        echo "The command failed"
        exit 1
    fi

    jj=`echo $output | grep 'Job'|tail -1|cut -f 2 -d "<"|cut -f 1 -d ">"`
    #check if it is a job array or not
    lines=`bjobs -r $jj|wc -l`
    if [[ "$lines" -gt 2 ]]; then
        echo "this seems an array job with $((lines-1)) jobs"

        vars=(`bjobs -A $jj |head -1`)
        vals=(`bjobs -A $jj |tail -1`)
        echo "$vars $vals"
        for ((i=0;i<${#vars[@]};i++)); do
            eval $(echo ${vars[$i]}=${vals[$i]})
        done
        if [ $PEND -gt 0 ]; then
            status="PEND"
        elif [ $RUN -gt 0 ]; then
            status="RUN"
        elif [ $DONE = $NJOBS ]; then
            status="DONE"
        elif [ $EXIT = $NJOBS ]; then
            status="EXIT"
        fi
    
        echo "will wait until job $jj finishes"
        while [ "$status" == "PEND" ] || [ "$status" == "RUN" ]; do
            echo "waiting job is $status"
            sleep 1
            vars=(`bjobs -A $jj |head -1`)
            vals=(`bjobs -A $jj |tail -1`)
            #assign vals to vars
            for ((i=0;i<${#vars[@]};i++)); do
                eval $(echo ${vars[$i]}=${vals[$i]})
            done
    
            if [ $PEND -gt 0 ]; then
                status="PEND"
            elif [ $RUN -gt 0 ]; then
                status="RUN"
            elif [ $DONE = $NJOBS ]; then
                status="DONE"
            elif [ $EXIT = $NJOBS ]; then
                status="EXIT"
            else
                status="Partial EXIT"
            fi
    
        done
        #sleep 130
        echo "job array $jj is finished with status $status"
        if [[ "$status" = "EXIT" ]]; then
            echo "will try to exit schedule"
            exit 1
        elif [[ "$status" = "Partial EXIT" ]]; then
            echo "It seems some of the jobs in the array have EXIT status"
        fi
    else
        echo "this seems a regular job with $((lines-1)) jobs"
        false
    fi
} || \
{ #try to see if it is a scheduled job on the cluster
    if [[ -z $check ]]; then
        echo "The command failed"
        exit 1
    fi
    jj=`echo $output | grep 'Job'|tail -1|cut -f 2 -d "<"|cut -f 1 -d ">"`

    status=`bjobs -r "$jj"|awk '{print $3}'|head -30|tail -1`
    echo "will wait until job $jj finishes"
    while [ "$status" == "PEND" ] || [ "$status" == "RUN" ]; do
            echo "waiting job is $status"
            sleep 1
            status=`bjobs -r "$jj"|awk '{print $3}'|tail -1`
    done
    sleep 130
    echo "singular job $jj is finished with status $status"
    if [ "$status" == "EXIT" ]; then
        echo "will try to exit schedule"
        exit 1
    fi
} || \
{ #Now check if the command failed 
#check=`echo $output | awk '{print $NF}'`
#echo "check is $check"
#check if $check is empty
if [[ -z $check ]]; then
    echo "The command failed"
    exit 1
fi
} || \
{ #catch seems either is is up to date or it is not running for some reason go to next step
echo "not sure what the command returned, not failed, not up to date and no job number ...check"
exit 1
}
} #end of function

function pipeline {
    if [ ! -f  16S_bacteria.cm ] || [ ! -f Makefile ] || [ ! -d scripts ] ; then
        echo "Hi, before starting the pipeline please make sure that you did the following:\n
        1. you copied the directory scripts/ into the current directory\n
        2. you copied the covariance file 16S_bacteria.cm to the current directory\n
        3. you copied Makefile to the current directory\n
        "
        echo "please do these configuration steps and try again"
        exit 1
    else
        continue
        echo "the pipeline is started"
        echo `date`
    fi
    local bigDB=$1
    local querySampleInfoF=$2
    local taxDB=$3
    local myemail=$4

    make logging
    make masked_reffasta.fasta bigDB="$bigDB"
    wait
    out_svs=`make querySVsDir querySampleInfoF="$querySampleInfoF"`
    checkJoboutput "$out_svs"

    if [ -f samplesMissed_SVs_multiplicity.csv ]; then
        nl=`cat samplesMissed_SVs_multiplicity.csv|wc -l`
        if [[ $nl -ge 2 ]]; then
            exit 1
        fi
    fi
    
    
    out_vs=`make vsearchOut/matched_*.fasta`
    checkJoboutput "$out_vs"
    #sleep 130 #seems that vsearch takes time to write results into fasta after finishing
    
    out_cvs=`make db_vsearch.fasta`
    checkJoboutput "$out_cvs"
    
    
    out_ncbi=`make dbvs.fasta`
    checkJoboutput "$out_ncbi"
    
    out_refal=`make dbvs.clean.align.fasta`
    checkJoboutput "$out_refal"
    
    out_rax=`make dbvs.clean.align.rax.fasta`
    checkJoboutput "$out_rax"
    
    make alignedQuery/sample_*.db.qc.clean.align.fasta
    out_alncbi=`make dbvs.clean.align.ncbi.fasta`
    checkJoboutput "$out_alncbi"
    out_seqI=`make dbvs.SeqInfo.csv myemail="$myemail"`
    
    
    out_tree=`make RAxML_fastTreeSH_Support.conf.root.ref.tre`
    checkJoboutput "$out_tree"
    
    checkJoboutput "$out_seqI"
    make checkpoint_seqInfo 
    
    out_taxI=`make dbvs.TaxaInfo.csv taxDB="$taxDB"`
    checkJoboutput "$out_taxI"

    out_treeNWK=`make RAxML_fastTreeSH_Support.conf.root.ref.nwk`
    checkJoboutput "$out_treeNWK"

    out_refP=`make allsamples.refpkg`
    checkJoboutput "$out_refP"
    
    out_jpl=`make phyloPlace/sample_*.jplace`
    checkJoboutput "$out_jpl"
    out_chkjpl=`make checkpoint_jplaceNames`
    checkJoboutput "$out_chkjpl"
    make analysis/allsamples.alphaDiversity.csv querySampleInfoF="$querySampleInfoF"

    out_jpstats=`make analysis/*_pplaceStats.csv`
    checkJoboutput "$out_jpstats"
    
    out_db=`make analysis/*.db`
    checkJoboutput "$out_db"
    out_chkdb=`make checkpoint_dbnames`
    checkJoboutput "$out_chkdb"
    
    out_bestR=`make analysis/*_bestRankStats.csv`
    checkJoboutput "$out_bestR"
    
    make analysis/oneRankEachSV_keepBest.csv
    sleep 10
    out_mulTax=`make analysis/multipleRanksSameSVMul.csv`
    checkJoboutput "$out_mulTax"
    
    out_ranks=`make analysis/ranksForTaxids.csv taxDB="$taxDB"`
    checkJoboutput "$out_ranks"
    
    #this intentionally put at last as it seems to occasionally corrupt the first jplace file
    make analysis/allsamples.jplace
    
    make clean
    echo "the pipeline is finished"
    echo `date`
}
export LSB_JOB_REPORT_MAIL="N"
export -f checkJoboutput 
export -f pipeline

bsub -q "standard" -n 1 -P microbiome -J "schedule" -oo "logging/schedule.log" "pipeline '$bigDB' '$querySampleInfoF' '$taxDB' '$myemail'"
