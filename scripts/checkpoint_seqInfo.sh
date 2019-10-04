#!/bin/bash
#requires python
module load python/3.5.2
function func1 {
echo -e 'from __future__ import print_function\nimport pandas as pd\ndf=pd.read_csv("dbvs.SeqInfo.csv")\nif df[df.tax_id.isnull()].shape[0]>0: print("some tax_id values are null in the dbvs.SeqInfo.csv file",df[df.tax_id.isnull()])\nelse: print("checkpoint passed successfully");'|python
#print("some tax_id values are null in the dbvs.SeqInfo.csv file");"
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "checkseqInfo" -oo "logging/chkseqInfo.log" "func1"
