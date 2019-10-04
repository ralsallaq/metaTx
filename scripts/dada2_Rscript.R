#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
csvFile <- read.csv(args[1])
datapaths <- csvFile$Path
R1F<-csvFile$R1
R2F<-csvFile$R2
fnFs <- paste(datapaths,R1F, sep="/")
fnRs <- paste(datapaths,R2F, sep="/") 
dirNames <- csvFile$dirName 
#extract sample names 
sample.names <- csvFile$sampleName
print(sample.names)

sessionInfo()

library(dada2);packageVersion("dada2")
outpath<-"./"

filt_path <- file.path(outpath, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,trimLeft=c(20,20),truncLen=c(290,210), maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
print(head(out))
#"""----To track which reads were lost in the filtering step. That would have to be done by reading the raw and filtered files back in and comparing them.------"""

apply(out,2,mean) #calculate the mean of the columns
apply(out,2,min)
Frac.Out=out[,2]/out[,1]
mean(Frac.Out)
print(summary(Frac.Out))
#Learn the error rates before dereplication
errF <- learnErrors(filtFs, multithread=TRUE) #very slow set multithread=TRUE
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"errF.rds")
saveRDS(errR,"errR.rds")


#Dereplication: combine all identical sequences into "unique sequences" with corresponding "abundance"
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
saveRDS(derepFs,"derepFs.rds")
saveRDS(derepRs,"derepRs.rds")

##########Denoising######
###Sample Inference:We are now ready to apply the core sequence-variant inference algorithm to the dereplicated data.
#This method removes indel errors (denoising) but chimeras remain 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) #Infer the sequence variants in each sample
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#"""-----The output of this function contains a $map of each input dereplicated unique sequence to the index of its correponding denoised sequence variant e.g. dadaFs$StoolDCM0SCR$map----"""
saveRDS(dadaFs, "dadaFs.rds")
saveRDS(dadaRs, "dadaRs.rds")
#Merge paired reads: reducing spurious sequence variants further more!
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers,"mergers.rds")
seqtab<-makeSequenceTable(mergers)
saveRDS(seqtab,"seqtab.rds")
#inspect the distribution of seq length:
seqLenDist=table(nchar(getSequences(seqtab)))
saveRDS(seqLenDist,"seqLenDist.rds")

#remove chimeric sequences:
seqtab.nochim<-removeBimeraDenovo(seqtab, method='consensus', multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"seqtab_nochim.rds")
sum(seqtab.nochim)/sum(seqtab) #this is ~96%
write.csv(seqtab.nochim, "seqtab_nochim.csv")

####Track reads through the pipeline to check progress:
getN<-function(x) sum(getUniques(x))
#sapply applies getN to each row (sample) of the table
trackFs<-cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
trackRs<-cbind(out, sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
#the first two columns are from out (filtering)
colnames(trackFs)<-c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(trackFs)<-sample.names
colnames(trackRs)<-c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(trackRs)<-sample.names
saveRDS(trackFs,"trackFs.rds")
saveRDS(trackRs,"trackRs.rds")

######getting fasta file for SVs without re-replication :
require(ShortRead)
tosave_df<-data.frame(sample=0,SVs_Abundance=0,seqtab_colind=0)
tosave_df<-tosave_df[-c(1),] #a trick to get correctly formatted dataframe to rbind to
tosave_df_missed<-data.frame(sample=0,SVs_Abundance=0,seqtab_colind=0)
tosave_df_missed<-tosave_df_missed[-c(1),] #a trick to get correctly formatted dataframe to rbind to
dir.create(args[2])
for(i in seq(nrow(seqtab.nochim))) {
    #extract concise sample name
    #sname<-strsplit(rownames(seqtab.nochim)[i],"_")[[1]][1]
    sname<-rownames(seqtab.nochim)[i]
    multiplicity=rowSums(seqtab.nochim)[i] #sum of all abundances for that sample
    #SVs with abundance >0 for the sample
    ii<-which(seqtab.nochim[i,]>0) #named int
    #if there any SVs with multiplicity >0 save them
    if (length(ii) > 0) { 

    temp_df<-data.frame(sname, unname(seqtab.nochim[i,ii]),unname(ii))
    names(temp_df)<-c("sample","SVs_Abundance","seqtab_colind")
    tosave_df<-rbind(tosave_df, temp_df)
    #seqs <- names(seqtab.nochim[i,ii]) #only SVs with multiplicity>0
    seqs <- names(ii) #only SVs with multiplicity>0
    #create id 
    id <- paste0(sname," Multiplicity=",unname(seqtab.nochim[i,ii])," colind=",unname(ii))
    #output a .fasta file is the seqs are not null
    #if (! is.null(seqs)) {
    writeFasta(object = ShortRead(sread = DNAStringSet(seqs), 
               id = BStringSet(id)), file = paste0("./SVs_fasta/",rownames(seqtab.nochim)[i],".fasta"), width=20000)
    cat("done with", i,"\n")
    }
    temp_df2 <-data.frame(sname, unname(seqtab.nochim[i,ii]),unname(ii))
    names(temp_df2)<-c("sample","SVs_Abundance","seqtab_colind")
    tosave_df_missed<-rbind(tosave_df_missed, temp_df2)
               }
write.csv(tosave_df,"samples_SVs_multiplicity.csv")
write.csv(tosave_df_missed,"samplesMissed_SVs_multiplicity.csv")
print("successfully finished dada2 analysis in R: check the following:1. you have a fasta file for each sample (no missingnesses) 2. check the seqtab_nochim.csv that it exists and looks fine")
