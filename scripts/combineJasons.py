import argparse
from Bio import SeqIO, Entrez, pairwise2
from Bio.SeqRecord import SeqRecord
import re
import os, sys, glob
import pandas as pd
import numpy as np
import json


def main():
    """analyze pplacer output from sqlite3 file"""
    saveFile=sys.argv[1]
    files=sys.argv[2:]
    print(files)
    result=[]
    for f in files:
        with open(f, "r") as infile:
            try:
                jj = json.load(infile)
            except:
                print("error loading json file")
            result.append(jj)

    print("saving jason file")
    with open(saveFile,"w") as outfile:
        json.dump(result, outfile)

if __name__ == "__main__":
    main()
