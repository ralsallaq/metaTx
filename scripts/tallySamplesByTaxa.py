#!/usr/bin/env python 
""" 
This module demultiplex the sequence file generated by decard to separate fasta files speciesd by SAMPLENAME_R1_ID.fasta
"""
from __future__ import print_function
import argparse
from Bio import SeqIO, Entrez, pairwise2
from Bio.SeqRecord import SeqRecord
import re
import os, sys, glob
import pandas as pd
import numpy as np
import random
import uuid


def main():
    """analyze pplacer output from sqlite3 file"""
    saveFile=sys.argv[1]
    svSeqTabF=sys.argv[2]
    files=sys.argv[3:][0].split(" ")
    print(files)
    df_seqtab=pd.read_csv(svSeqTabF, index_col=0)
    #replace SVs by their colind as in R (start from 1 and end at the number of SVs)
    df_seqtab.columns = list(range(1,df_seqtab.shape[1]+1))
    print(df_seqtab.head())

    df_list=[]
    all_taxa=np.array([])
    all_taxa_tosave=np.array([])
    for f in files:
        sname=os.path.basename(f).split("_bestRankStats")[0]
        df1=pd.read_csv(f,index_col=0)
        df1.loc[:,'tax_id']=df1.loc[:,'tax_id'].astype(str)
        df1.loc[:,"Multiplicity"] = df1.name.apply(lambda row: float(row.split("_")[-3])) #the third from the last
        df1.loc[:,"colind"]=df1.name.apply(lambda row: int(row.split("_")[-1])) #this is the distinctive index of SVs in seqtab.nochim (not pythonic index), think of it as an id for SVs
        unique_SVs=df1["colind"].unique()
        ind_unique=~df1["colind"].duplicated()
        totalMultiplicitySample=df1.loc[ind_unique,"Multiplicity"].sum()
        print(totalMultiplicitySample)
        #group by distinct SVs
        gpby=df1.groupby('colind', as_index=False)
        temp = gpby.first()[['rank','tax_id','tax_name','Multiplicity', 'colind','likelihood']]
        assert(temp['Multiplicity'].sum()==totalMultiplicitySample), "the multiplicities of SVs are not unique in temp data frame" 
        assert(temp['Multiplicity'].sum()==df_seqtab.loc[sname,:].sum()), "the total multiplicity for the sample file "+f+" is not equal to that in the seqtab!"
        #print(list(gpby)[2][1])
        #print(list(gpby)[3][1])

        """ if for each group distinct SVs exist (distinct colind values) then the following would work
        however this is not the case for every taxonomy
        temp.set_index('tax_name',inplace=True)
        multiplicity=gpby.sum()[['tax_name','Multiplicity']]
        multiplicity.set_index('tax_name',inplace=True)
        temp = temp.merge(multiplicity, how="outer", left_index=True, right_index=True)
        temp = temp.reset_index()
        temp.rename(columns={'Multiplicity':sname},inplace=True) """

        #each SVs can be assigned different placements we pick the one with the highest likelihood
        for rowi,row in temp.iterrows():
            gp = list(gpby)[rowi][0]
            #print(gp, row.values) these are equal but with different order
            df_gp = list(gpby)[rowi][1]
            #save all taxids in array
            all_taxa_tosave = np.append(all_taxa_tosave, df_gp['tax_id'].astype(str).values)
            #more than one different taxonomies for the same SV
            if df_gp.shape[0]>1: #choose the highest likelihood placement:
                #print("before",temp.loc[rowi,:])
                ind=df_gp['likelihood'].idxmax()
                taxids_sameSV = ",".join(df_gp['tax_id'].astype(str).values)


                tax_names_sameSV_ = df_gp.loc[:,'tax_name'].apply(lambda row: str(row).replace("]","").replace("[","").split(" "))
                tax_names_sameSV = set(tax_names_sameSV_.iloc[0])
                for i in range(1,tax_names_sameSV_.shape[0]):
                    tax_names_sameSV = tax_names_sameSV.union(set(tax_names_sameSV_.iloc[i]))
               

                #first assign the default
                temp.loc[rowi,:]=df_gp.loc[ind,['rank','tax_id','tax_name','Multiplicity', 'colind','likelihood','sname']]
                #then change to add multiple ranks for each SV
                temp.loc[rowi,'tax_id']=taxids_sameSV
                temp.loc[rowi,"tax_name"]="/".join(sorted(tax_names_sameSV)) #sort to get the capitalized word first
                #temp.loc[rowi,"tax_name"]="/".join(tax_names_sameSV) #sort to get the capitalized word first
                #print("after",ind,temp.loc[rowi,:])
            else:
                #remove square brakets from names 
                temp.loc[rowi,'tax_name']=df_gp['tax_name'].astype(str).values[0].replace("]","").replace("[","")
                continue

        #another grouping is necessary to sum multiplicities under the same taxonomies for different SVs
        gpbyTaxa=temp.groupby(['rank','tax_id','tax_name'], as_index=False)
        temp=gpbyTaxa.first()
        temp.loc[:,'Multiplicity']=gpbyTaxa.sum()['Multiplicity']

        temp.rename(columns={'Multiplicity':sname},inplace=True)
        temp.drop(['colind','likelihood'],axis=1, inplace=True)
        all_taxa = np.append(all_taxa, temp['tax_id'].astype(str))
        df_list.append(temp)

        print(temp[sname].sum())
        print("-------------------")

    all_taxa=np.unique(all_taxa) 
    print("number of distinct taxa across all samples is:", all_taxa.shape[0])
    np.save("analysis/allTaxIds.npy",all_taxa_tosave)
    
    #combine the samples and tally them against taxa
    #print(df_list[1].sum()[sname])
    df_all=df_list[0]
    for df in df_list[1:]:
        df_all = df_all.merge(df, how="outer", on=['rank','tax_id','tax_name'])

    df_all.to_csv(saveFile) #sum of the columns for each sample would give the total multiplicity of SVs in the seqtab.nochim dada2 table i.e. sum of column for sample1 = sum(seqtab.nochim[sample1,])

if __name__ == "__main__":
    main()

