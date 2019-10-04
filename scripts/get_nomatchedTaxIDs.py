#!/usr/bin/env python 
""" 
This module compares the unique taxids from the aligned fasta file  and the unique taxids from the taxonomyInfo file and output a list of no matches
These should be inquired on https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi and be replaced in the aligned ref and aligned query/ref files
to match primary taxid (normally found in the taxonomyInfo file)
"""
from __future__ import print_function
import argparse
from Bio import SeqIO, Entrez, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import re
import os, sys, time
import pandas as pd
import numpy as np
from Bio import Entrez
import tempfile
from shutil import copy2
######################################################################


p_desc=re.compile(r"""
(description)   #the keys group
(=")  #the values starter
([^"]+) #the values: everything without the = sign even one time
(")  #the values terminator
""", re.VERBOSE)


p_taxa=re.compile(r"""
(taxonomy)   #the keys group
(=")  #the values starter
([^"]+) #the values: everything without the = sign even one time
(")  #the values terminator
""", re.VERBOSE)

p_organism=re.compile(r"""
(organism)   #the keys group
(=")  #the values starter
([^"]+) #the values
(")  #the values terminator
""", re.VERBOSE)

p_taxid=re.compile("""
        (ncbi_tax_id)
        (=")
        (.*)
        (")
""", re.VERBOSE)


#######Functions to fetch taxonomy from NCBI db #############
Entrez.email = "ralsalla@stjude.org"
if not Entrez.email:
    print ("you must add your email address")
    sys.exit(2)

def Entrez_NCBI():
    """ This function displays all searchable fields available in Entrez"""
    data = Entrez.read(Entrez.einfo(db="taxonomy"))
    for field in data["DbInfo"]["FieldList"]:
        print("%(Name)s, %(FullName)s,%(Description)s" % field)

def get_taxid(species):
    """to get data from ncbi taxomomy, we need to have the taxid.  we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList'][0]

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)

def get_lineage(data):
    """once you have the data from get_tax_data fetch the lineage"""
    uprank=['kingdom','phylum','class','order','family','genus']
    lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in uprank}
    try: 
        lineage['kingdom']
    except KeyError:
        dict_ = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank']=="superkingdom"}
        dict_['kingdom']=dict_.pop('superkingdom')
        lineage.update(dict_)
    return lineage



def main():
    """converts fasta from alignment clean titles to fasta with standard titles in ncbi format"""

    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--input_fasta','-i1', help='aligned ref or ref/query fasta file with clean format titles (output from aligners)', required=True)
    args_parser.add_argument('--input_primarytaxa','-i2', help='csv taxonomy information file with primary taxids', required=True)
    args_parser.add_argument('--out_nomatch','-o', help='csv file with heads tax_id, primary_tax_id species_name', required=True)
    args = args_parser.parse_args()

    # First handle the files
    input_fasta = args.input_fasta
    input_taxa = pd.read_csv(args.input_primarytaxa)
    output_file = args.out_nomatch

    with open(input_fasta,'rt') as subrecs:
        tax_id=[]
        species_name=[]
        df_alignFTaxa=pd.DataFrame()
        for title, sequence in SimpleFastaParser(subrecs):
            ttl=title.split(None,1)[0]
            desc=title.split(None,1)[1]
            try:
                taxid=p_taxid.findall(desc)[0][2]
                organism= p_organism.findall(desc)[0][2].replace("=","").replace("-","_").replace("[","").replace("]","").replace(".","_")
                #print(taxid, organism)
            except:
                taxid=""
                name = ttl.split(";")[0]
                print(name, desc)
                Len =  ttl.split(";")[1].split("=")[1]
                organism=""
                taxonomy=""
                details=""
            tax_id.append(taxid)
            species_name.append(organism)

        df_alignFTaxa.loc[:,'tax_id']=tax_id
        df_alignFTaxa.loc[:,'tax_id'] = df_alignFTaxa.loc[:,'tax_id'].astype(int)
        df_alignFTaxa.loc[:,'species_name']=species_name
        df_alignFTaxa.drop_duplicates(inplace=True)
        df_nomatch=df_alignFTaxa[~df_alignFTaxa.tax_id.isin(input_taxa.tax_id.astype(int).values)]
        print(df_alignFTaxa[~df_alignFTaxa.tax_id.isin(input_taxa.tax_id.astype(int).values)])
        primary_taxid=[]
        for rowi, row in df_nomatch.iterrows():
            recordFitched = False
            while not recordFitched:
                try:
                    data=get_tax_data(str(row.tax_id))
                    print(data[0]['TaxId'])
                    primary_taxid.append(data[0]['TaxId'])
                    recordFitched = True
                except:
                    time.sleep(8)
                    recordFitched = False

        df_nomatch.loc[:,'primary_tax_id']=primary_taxid
        df_nomatch.to_csv(output_file, index=False)

    print("csv file successfully saved")


    #savefasta(df_merge, no_qc_outf)

if __name__ == "__main__":
    main()
