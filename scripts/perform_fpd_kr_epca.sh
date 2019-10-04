module load python/3.5.2
module load pplacer
#This bash script executes the phylogenetic placement using pplacer  
if [ ! -d "phyloPlace" ]; then
    echo "seems you have not done the pplacer step"
    exit 1
fi
if [ ! -d "analysis" ]; then
    mkdir analysis
fi

refpkgdir="$1"
querySampleInfoF="$2"
shift
shift
#the rest are all the files including the one to be excluded 
allfiles=$@
#get files that we will include in epca analysis
include1=`python -c "from __future__ import print_function;import pandas as pd; df=pd.read_csv('$querySampleInfoF', index_col=0); print(df[df['Eval']==1]['sampleName'].values)"`
xclude1=`python -c "from __future__ import print_function; import pandas as pd; df=pd.read_csv('$querySampleInfoF', index_col=0); print(df[df['Eval']==0]['sampleName'].values)"`
snamesI=`echo $include1|cut -d "[" -f2 | cut -d "]" -f1`
snamesX=`echo $xclude1|cut -d "[" -f2 | cut -d "]" -f1`
#number of sample files
nAll=`echo $allfiles|wc -w`
nFs=`echo $snamesI|wc -w`
NToX1=$((nAll-nFs))
NToX2=`echo $snamesX|wc -w`

if [[ "$NToX1" -ne "$NToX2" ]]; then

    echo "$NToX1 is not equal to $NToX2 ! check the files to be excluded, it might be that some samples did not have SVs in dada2 step, if this is the case remove these samples from the file $querySampleInfoF"
    exit 1
fi

files=() #define an empty array
for ((i=1;i<="$nFs";i++)); do
    var=`echo $snamesI | cut -f $i -d " "|cut -f 2 -d "'"`
    echo $var
    #append to the array
    files+=("phyloPlace/sample_"$var".jplace")
done
#number of sample files
echo "ref packag=" $refpkgdir 
echo "number to x=" $NToX1 
echo "samples to x=" $snamesX 
echo "files for the analysis are the files rendered by the wildcard from Make as make will automatically exclude explicitly specified files from the wildcard" ${files[@]}
echo "number of files to include in analysis=" $nFs
echo "number of all files=" $((nFs+NToX1))
#test the files
check=`ls -p ${files[@]}|wc -w`
if [[ "$check" -ne "$nFs" ]]; then
    echo "check the jplace files ${files[@]} seems not all of them exist"
    exit 1
fi

analysisFs=`ls -p ${files[@]}`

##go through the sample files
function func1 {

    files=$1
    refpkgdir=$2

    echo $refpkgdir $files 
    echo "-----------"
    echo "get alpha diversity in terms of Jost Lou effective numbers of order 0 (SV richness), 1 (Shannon effective number), and 2 (Simpson's effective number)"
    guppy fpd --pp --csv --chao-d 0,1,2 --theta 0,1 -o analysis/allsamples.alphaDiversity.csv $files

    echo "calculate the generalized weighted UniFrac distance (Kantorovich-Rubinstein distance) between samples using the collection of placements in each sample's jplace file on the ref tree" 
    guppy kr --point-mass --pp --list-out -c $refpkgdir -o allsamples.unifrac.csv --out-dir analysis/ --seed 1 $files 

    echo "perform edge principle components on the files"
    guppy epca --point-mass --pp --som 2 --kappa 1 --prefix epcaOUT --out-dir analysis/ -c $refpkgdir $files

}

export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 4 -R "span[hosts=1] rusage[mem=6000]"  -P microbiome -J "fpd_kr_epca" -oo "logging/perform_epca.log" "func1  '$analysisFs' '$refpkgdir'"
