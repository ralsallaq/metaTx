#!/usr/bin/env /hpcf/apps/python/install/3.5.2/bin/python
import sys, glob
from os import listdir
from os.path import isdir, isfile, join
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Entrez
import numpy as np
import pandas as pd
import sqlite3
import argparse
import json #to read jplace files
from ast import literal_eval


Entrez.email = "ralsalla@stjude.org"

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)

def get_lineage(data):
    """once you have the data from get_tax_data fetch the lineage"""
    uprank=['kingdom','phylum','class','order','family','genus','species']
    lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in uprank}
    try:
        lineage['kingdom']
    except KeyError:
        try:
            dict_ = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank']=="superkingdom"}
            dict_['kingdom']=dict_.pop('superkingdom')
            lineage.update(dict_)
        except KeyError:
            lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx']}
    return lineage


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except:
        print("Error! connection could not be established")
        sys.exit(1)

    return None


def get_lineage_ids(taxid, conn):
    """ This function gets the names and ids if all parents of the given id """ 

    query="SELECT nd.tax_id, nd.parent_id, nd.rank, na.tax_id, na.tax_name, na.name_class from nodes nd inner join names na on nd.tax_id=na.tax_id where na.name_class=='scientific name' AND na.tax_id==" + "'"+taxid+"'"

    df = pd.read_sql_query(query, conn)
    #print(df)

    df.columns=['tax_id', 'parent_id', 'rank', 'tax_id_drop', 'tax_name', 'name_class']
    df.drop("tax_id_drop",axis=1, inplace=True)
    rankorder=np.array(['no_rank','superkingdom','phylum','class','order','family','genus','species'])[::-1]
    #print(df['rank'])
    if not df['rank'].iloc[0].strip() in rankorder:
        rankorder=np.append(rankorder,df['rank'].iloc[0])
    rank_ind=np.where(df['rank'].iloc[0]==rankorder)[0][0]
#    if len(df.tax_name.iloc[0].split(" "))>=2:
#        lineage={rankorder[rank_ind]:" ".join(df.tax_name.iloc[0].split(" ")[1:])}
#    else:
#        lineage={rankorder[rank_ind]:df.tax_name.iloc[0]}

    lineage={rankorder[rank_ind]:df.tax_name.iloc[0]}
    ids={rankorder[rank_ind]: taxid}
    stop=False
    temp=df.copy()
    #print(lineage, ids)
    while not stop:
        parent_id=temp['parent_id'].iloc[0]
        if parent_id is None or parent_id=="" or parent_id=='0':
            stop=True
            break
        #print(parent_id)
        query="SELECT nd.tax_id, nd.parent_id, nd.rank, na.tax_id, na.tax_name, na.name_class from nodes nd inner join names na on nd.tax_id=na.tax_id where na.name_class=='scientific name' AND na.tax_id==" + "'"+parent_id+"'"
        temp = pd.read_sql_query(query, conn)
        temp.columns=['tax_id', 'parent_id', 'rank', 'tax_id_drop', 'tax_name', 'name_class']
        temp.drop("tax_id_drop",axis=1, inplace=True)
        lineage.update({temp['rank'].iloc[0]:temp['tax_name'].iloc[0]})

        ids.update({temp['rank'].iloc[0]: temp.tax_id.iloc[0]})

        

    return lineage, ids


def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--sqlitedb','-db', help='sqlite3 file with taxonomy data', required=True)
    args_parser.add_argument('--taxTable','-l', help='csv file with a column name "tax_id" encompassing a list of taxids to fetch their taxonomy', required=True)
    args_parser.add_argument('--output','-o', help='the taxids and their fetched taxonomies output file', required=True)
    args = args_parser.parse_args()

    sqlF = args.sqlitedb
    conn = create_connection(sqlF)
    taxids_notsplit=pd.read_csv(args.taxTable)['tax_id'].astype(str)
    taxids=set()
    for rowi, row in taxids_notsplit.iteritems():
        taxids=taxids.union(set(row.split(",")))
    taxids=list(taxids)

    print(taxids_notsplit.shape[0], len(taxids))
    outputF=args.output
    #ranks of interest
    ranks=np.array(['superkingdom','phylum','class','order','family','genus','species'])
    with conn:
        #get the names of tables:
#        query="SELECT name FROM sqlite_master WHERE type='table';"
#        table_names = pd.read_sql_query(query, conn)
#        # ranks, source, nodes, names, merged
#        print("table names are:\n", table_names, "\n")
#        query="SELECT * from names"
#        print(pd.read_sql_query(query, conn).head())
#        query="SELECT * from nodes"
#        print(pd.read_sql_query(query, conn).head())

        taxa=pd.DataFrame({'lineage':np.nan,'ids':np.nan},index=taxids)
        tofetch=[]
        for rowi, row in taxa.iterrows():
            print(rowi)

            try:
                lineage, ids=get_lineage_ids(rowi, conn)
                taxonomy=[]
                for r in ranks:
                    try:
                        taxonomy.append(lineage[r])
                    except:
                        pass
                #taxonomy=np.array([lineage[k] for k in lineage.keys() if k in ranks])[::-1]
                taxa.loc[rowi,'lineage']=[lineage]
                #print(taxonomy)
                taxa.loc[rowi,'taxonomy']=[";".join(taxonomy)]
                taxa.loc[rowi,'ids']=[ids]
            except:
                #tofetch.append(str(rowi))
                lineage = get_lineage(get_tax_data(str(rowi)))
                taxonomy=[]
                for r in ranks:
                    try:
                        taxonomy.append(lineage[r])
                    except:
                        pass
                #taxonomy=np.array([lineage[k] for k in lineage.keys() if k in ranks])[::-1]
                taxa.loc[rowi,'lineage']=[lineage]
                taxa.loc[rowi,'taxonomy']=[";".join(taxonomy)]
                taxa.loc[rowi,'ids']=[rowi]
                #print("except",[lineage[r] for r in ranks]) 

#    tofetch="/n".join(ids)
#    ids="/n".join(ids)
#    lineages = get_lineage(get_tax_data(str(rowi))) 
#    if len(tofetch)>0:
#        print(lineages, "\n",ids)
    print(taxa)
    #taxa['taxonomy'].apply(lambda row: row[0]).to_csv(outputF)
    taxa['taxonomy'].apply(lambda row: row[0] if type(row) is list else row).to_csv(outputF)
    return taxa
            

if __name__ == '__main__':
        main()

