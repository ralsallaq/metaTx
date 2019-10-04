#!/usr/bin/env 
from __future__ import print_function 
import sys, glob
from os import listdir
from os.path import isdir, isfile, join
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
import sqlite3
import argparse
import json #to read jplace files


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
                                         
    return None

def select_lineage_for_queries(conn, desired_rank='species'):
        """
        Query all rows in the tasks table
        :param conn: the Connection object
        :return:
        """
        cur = conn.cursor()
        querry_values = (desired_rank)
        cur.execute("select placement_id, tax_id, likelihood FROM multiclass as mc \
                     JOIN placement_classifications AS pc USING (placement_id) \
                     JOIN taxa USING (tax_id) \
                     JOIN ranks USING (rank) \
                     WHERE pc.rank = ? AND pn.name = ?",querry_values)
                 
        rows = cur.fetchall()
                   
        for row in rows:
            print(row)


def main():
    """analyze pplacer output from sqlite3 file"""

    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--sqlitedb','-db', help='sqlite3 file with classification data', required=True)
    args_parser.add_argument('--queryFasta','-q', help='fasta file for query', required=True)
    args_parser.add_argument('--outputCSV','-oc', help='CSV file with classification data for each query', required=True)
    args = args_parser.parse_args()
    sqlF = args.sqlitedb
    queryF = args.queryFasta
    outf = args.outputCSV
    desired_rank='species'
    with open(queryF,'rt') as q_hdl:
        q_ttls=[]
        for title, sequence in SimpleFastaParser(q_hdl):
            full_ttl=title
            q_ttls.append(title.split(None,1)[0].strip())

    with open(outf,"wt") as out_hdl:
        conn = create_connection(sqlF)
        with conn: #use context manager so you do not have to commit
            #get the name of the tables:
            query="SELECT name FROM sqlite_master WHERE type='table';"
            table_names = pd.read_sql_query(query, conn)
            print("table names are:\n", table_names, "\n")
            query="SELECT mc.name, mc.rank, mc.tax_id, taxa.tax_name,mc.likelihood from multiclass as mc join taxa on taxa.tax_id==mc.tax_id where want_rank=="+"'"+desired_rank+"'"
            query_classifications = pd.read_sql_query(query,conn)
            print("classifications for query are \n", query_classifications.head())
        query_classifications = query_classifications[query_classifications.name.isin(q_ttls)]
        #check if the desired rank is achieved and how many ranks off 
        query_classifications.loc[:,'achieved_rank'] = query_classifications['rank'].apply(lambda r:r.split("_")[0])
        #ranks off : 0 (desired rank achieved), 1 (one rank off),..etc
        query_classifications.loc[:,'ranks_off']=np.where(query_classifications['achieved_rank'] == desired_rank,0, np.where(query_classifications['achieved_rank']=='genus',1, \
                np.where(query_classifications['achieved_rank']=='family',2,np.where(query_classifications['achieved_rank']=='order',3,np.where(query_classifications['achieved_rank']=='class',4, \
                np.where(query_classifications['achieved_rank']=='phylum',5,np.where(query_classifications['achieved_rank']=='superkingdom',6,np.nan)))))))
        print(query_classifications.head())
        query_classifications.to_csv(out_hdl, index_label='index')


if __name__ == '__main__':
    main()

