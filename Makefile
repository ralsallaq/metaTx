SHELL :=/bin/bash

#Dummy targets:
all: logging masked_reffasta.fasta ./vsearchOut/matched_*.fasta db_vsearch.fasta dbvs.fasta dbvs.clean.align.fasta dbvs.clean.align.rax.fasta RAxML_fastTreeSH_Support.conf.root.ref.tre  alignedQuery/sample_*.db.qc.clean.align.fasta allsamples.refpkg analysis/multipleRanksSameSV.csv 

clean: 
	rm --force temp_* temp1_*

.PHONY: all clean

# """" create logging directory for log files
logging:
	mkdir logging

# """" link bigDB 
masked_reffasta.fasta: 
	rm masked_reffasta.fasta; ln -s $(bigDB) $@

# """" run dada2 on raw sequences
querySVsDir svSeqTab: ./scripts/dada2_Rscript.R $(querySampleInfoF) 
	./scripts/get_dada2_svs.sh  $^ $(PWD)/SVs_fasta
	wait
	ln -s $(PWD)/SVs_fasta querySVsDir; ln -s $(PWD)/seqtab_nochim.csv svSeqTab

# """" vsearch masked_ref for 80% similarity with query reads
vsearchOut/matched_%.fasta: masked_reffasta.fasta querySVsDir/%.fasta
	./scripts/vsearch_reffasta_forQueries.sh $^

# """" combine and dereplicate vsearch results 
db_vsearch.fasta: 
	./scripts/combine_derep_vsearchFasta.sh 

# """" add taxonomy to combined dbvsearch, note dbvs.fasta has now new ncbi taxids:
dbvs.fasta: db_vsearch.fasta masked_reffasta.fasta 
	./scripts/assignTaxonomyToVsearchFasta.sh $^ $@

# """" clean and align ref fasta:
dbvs.clean.align.fasta: dbvs.fasta 
	./scripts/cleanAlignRef.sh $^


# """" prep alignment file for rax
dbvs.clean.align.rax.fasta: dbvs.clean.align.fasta 
	./scripts/prepForRaxTreeBuild.sh $<
# build a reference tree: 
RAxML_fastTreeSH_Support.conf.root.ref.tre: dbvs.clean.align.rax.fasta 
	./scripts/buildRefTree.sh  $^ 

    
# """" fix the titles in the alignment ref fasta
dbvs.clean.align.ncbi.fasta raxml_heads ncbi_heads: dbvs.clean.align.rax.fasta dbvs.fasta 
	./scripts/fromCleanToStandardTitle.sh $^

# """" get SeqInfo and fasta files that only has primary taxids (no old taxids)
dbvs.SeqInfo.csv: dbvs.clean.align.ncbi.fasta 
	./scripts/get_SeqInfoCSV.sh $^ $(myemail) $@

# """" checkpoint  make sure that seqInfo do not have Null taxids
checkpoint_seqInfo:
	./scripts/checkpoint_seqInfo.sh

# """" save a copy of the primarytaxid fasta:
dbvs.clean.align.ncbi.primaryTaxids.fasta_save:
	cp dbvs.clean.align.ncbi.primaryTaxids.fasta dbvs.clean.align.ncbi.primaryTaxids.fasta_save

# """" get TaxInfo file
dbvs.TaxaInfo.csv: dbvs.SeqInfo.csv $(taxDB)
	./scripts/get_TaxaInfoCSV.sh $^ $@


# """" change tip labels in the ref tree to sequence ids as in the seqInfo file 
RAxML_fastTreeSH_Support.conf.root.ref.nwk: RAxML_fastTreeSH_Support.conf.root.ref.tre raxml_heads ncbi_heads 
	./scripts/fromCleanToSeqID_RAxMLTree.sh $^ 

# """" create refpkg

allsamples.refpkg: dbvs.SeqInfo.csv dbvs.TaxaInfo.csv RAxML_fastTreeSH_Support.conf.root.ref.nwk RAxML_info.ref.tre dbvs.clean.align.ncbi.primaryTaxids.fasta 
	./scripts/createRefpkgWtables.sh $^ 


# """" align query files : to make it just run make alignedQuery/sample_*.db.qc.clean.align.fasta
alignedQuery/sample_%.db.qc.clean.align.fasta: dbvs.clean.align.rax.fasta  querySVsDir/%.fasta 
	./scripts/cleanAlignQuery_bysample.sh $^ 

# """" place query reads on the reference tree (note that this step will take care of correcting titles of the reference sequences in the query alignments to NCBI titles) : to make it just run make phyloPlace/sample_*.jplace""""
phyloPlace/sample_%.jplace: allsamples.refpkg alignedQuery/sample_%.db.qc.clean.align.fasta 
	./scripts/phyloPlace.sh $^

#checkpoint make sure that the sample names are the same in the jplace and in the alignment files so that next steps are good to go
checkpoint_jplaceNames:
	./scripts/checkpoint_jplace.sh

# """" analyze placements using guppy run by: make analysis/*.db; 
analysis/%.db analysis/%.csv analysis/%.phyloxml analysis/%.edpl.csv  analysis/%.adcl.csv: allsamples.refpkg alignedQuery/sample_%.db.qc.clean.align.ncbi.fasta phyloPlace/sample_%.jplace 
	./scripts/prep_tovisualize.sh $^

#checkpoint 3 check that there is a correspondence between sequence Ids in *.db files and the *.assembled.qc.clean.fasta files
checkpoint_dbnames:
	./scripts/checkpoint_dbnames.sh

# """" get best rank stats
analysis/%_bestRankStats.csv: analysis/%.db alignedQuery/%.assembled.qc.clean.fasta
	./scripts/bestRankForReads.sh $^

# """" tally samples by taxa
analysis/multipleRanksSameSVMul.csv: svSeqTab analysis/*_bestRankStats.csv 
	./scripts/tallySamplesByTaxa.sh $@ $^
# """" create a file for the ranks for each taxid in the taxTable
analysis/ranksForTaxids.csv: analysis/multipleRanksSameSV.csv $(taxDB) 
	./scripts/get_taxa_fromTaxaDB.sh $@ $^
	
# """" tally samples by taxa to create a taxTable but only include highest likelihood
analysis/oneRankEachSV_keepBest.csv: svSeqTab analysis/*_bestRankStats.csv
	./scripts/tallySamplesByTaxa_keepBest.sh $@ $^

# """" Get placement richness and distal length as pplacer stats
analysis/%_pplaceStats.csv: phyloPlace/sample_%.jplace
	./scripts/jplaceStats.sh $^ $@
	
# """" alpha diversity, beta diversity, Edge PCA and splitify matrix:
analysis/allsamples.alphaDiversity.csv analysis/allsamples.unifrac.csv analysis/epcaOUT.trans analysis/epcaOUT.proj analysis/epcaOUT.xml analysis/epcaOUT.som analysis/epcaOUT.som.xml: allsamples.refpkg $(querySampleInfoF) phyloPlace/sample_*.jplace 
	./scripts/perform_fpd_kr_epca.sh $^

#"""" combine jplace files into one:
analysis/allsamples.jplace: phyloPlace/sample_*.jplace
	./scripts/combineJasons.sh $@ $^
