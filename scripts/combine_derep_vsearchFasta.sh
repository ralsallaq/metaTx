#!/bin/bash
#requires seqmagick
module load python/2.7.13
files=`ls -p vsearchOut/matched_*.fasta`
#number of sample files
nFs=`echo $files |wc -w`
echo "will combine and dereplicate $nFs files"
echo "cleaning temporary files that are used by the script"
rm -f tempF temp2.fasta

function func1 {
local file=$1
local tempF=$2
    echo "------processing : $file-----------------"
    cat $tempF $file > temp2.fasta
    #deduplicate sequences inplace
    wait
    seqmagick mogrify --deduplicate-sequences temp2.fasta
    wait
    mv -f temp2.fasta $tempF
    echo "------done with  $file --------"
    echo "$tempF"

}

function func2 {
    local files=$1
##python way
cat <<EOF > ./scripts/temp.py
#!/usr/bin/python
from __future__ import print_function
import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
outputF=sys.argv[1]
files=sys.argv[2:][0].split("\n")
print(files)
with open(outputF,'wt') as out_hdl:
    ttls=[]
    seqs=[]
    for f in files:
        print("processing file ", f)
        with open(f, 'rt') as inp_hdl:
            for title, sequence in SimpleFastaParser(inp_hdl):
                if (title in ttls) and (sequence in seqs):
                   continue
                ttls.append(title)  
                seqs.append(sequence)
    for ttl, seq in zip(ttls, seqs):
        out_hdl.write(">%s\n%s\n" % (ttl, seq))
EOF
python ./scripts/temp.py db_vsearch.fasta "$files"

 echo "cleaning temporary files that are used by the script"
 rm -f ./scripts/temp.py
}

export LSB_JOB_REPORT_MAIL="N"
export -f func1 
export -f func2

bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=5000]" -P microbiome -J "combineVFs" -oo "logging/combine_derep_vsearchFasta.log" "func2 '$files'"
