Acknowledgement:
This pipeline was not possible without the help of two individuals:
1. Elisa Margolis
2. Jonathan Golob who also has a similar pipeline that uses dada2+pplacer but utilizes  nextflow instead of Gnu Make. In fact these two pipelines were developed in tandem and also their outputs were compared.  

=====================================================================================
    
metaTx is an automated pipeline that is designed for now to run on high performance computational clusters (HPC) using LSF (bsub/IBM)
It uses numerous shell scripts that are designed to achieve most of the tasks in parallel on the HPC. The flow dependencies are controlled by using gnu make and automated scheduling by using a shell script that is designed specifically for make/lsf. The shell scripts including the scheduling one can be easily (hopefully) redesigned for any batch system such as slurm or torque.
Currently the pipeline is optimized to run for 2x300 alumina pair-end reads with a primer of 20 bps or so. For other cases the script scripts/dada2_Rscript.R should be edited according to dada2 recommended settings for the specific case. 

Requirements:
To run the pipeline you need to have a couple of programs available in your path on the HPC and these are listed per each shell script in the requirements.txt file 
These programs might be available on the HPC cluster (check that by module avail PROGRAMNAME) or can be installed into conda virtual environments under user account 
At the top of the scripts if the module available in the HPC it will be called into the script that needs it via:
module load PROGRAMNAME/PROGRAMVersion
If it is available in user conda env then it will be called via: 
activate=path/to/bin/activate
source $activate PROGRAMNAME
The scripts for metaTx use the module-load track as many of the programs were installed into the HPC at St.Jude  

===================================================================================

How to run:
To run the pipeline (1)  prepare a csv file with each row corresponding to a distinctly named sample with the following required fields:
Path=full path for R1/R2 reads (assuming both reads for the sample are under this directory)  
R1=full name of R1 file for the sample (without the path)
R2=full name of R2 file for the sample (without the path)
sampleName=the desired distinct name for the sample this name will be kept for downstream analysis and it must be distinct for each sample
Eval=a binary indicator to indicate whether to keep (value =1) the sample for edge principle component analysis or to exclude (value=0)
Other metadata columns of interest 
(2) copy two files and a directory into the directory in which you wish to run the pipeline: 16S_bacteria.cm  Makefile and scripts/ respectively
(3) run the pipeline: 
./scripts/schedule.sh -r /absolute/path/to/HugeRefFasta.fasta -f /absolute/path/to/sampleInfo.csv -d /absolute/path/to/taxonomy.db
were:
/absolute/path/to/HugeRefFasta.fasta: path to the reference fasta file of the database to be used (RDP, NCBI)
/absolute/path/to/sampleInfo.csv: the sample info csv file previously prepared  
/absolute/path/to/taxonomy.db: the location of the NCBI taxonomy sqlite3 database containing taxonomy information on the taxids in the  HugeRefFasta.fasta fasta file
(4) if something goes wrong check the schedule.log file under the logging directory created by the pipeline and the log file of the last step at which the pipeline breaks
(5) In case of a batch system different than LSF (e.q. Moab, Slurm, Torque) then one needs to change the shell scripts including the schedule.sh accordingly. please post in the issues if you need help with this.    
