#!/usr/bin/env python 
""" 
This module changes the tips on the raxml tree to seqIDs
"""
from __future__ import print_function
import argparse
from Bio import SeqIO, Entrez, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import re
import os, sys
import pandas as pd
import numpy as np
from Bio import Entrez
import tempfile
from shutil import copy2
######################################################################
def main():
    """replaces tip labels on raxml newick tree from raxml format to SeqIDs as in the SeqInfo file"""

    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--input_raxmlH','-i1', help='text file with clean format titles (obtained via grep "^>" ncbi_ncbi.clean.align.taxa.rax.fasta|cut -f2 -d ">")', required=True)
    args_parser.add_argument('--input_ncbiH','-i2', help='text file with ncbi format titles (obtained via grep "^>" ncbi_ncbi.taxa.fasta|cut -f2 -d ">")', required=True)
    args_parser.add_argument('--input_tree','-iT', help='raxml tree file with raxml acceptable tips', required=True)
    args_parser.add_argument('--out_tree','-oT', help='newick tree file with standard ncbi seqIDs with . replaced by _', required=True)
    args = args_parser.parse_args()

    # First handle the files
    with open(args.input_raxmlH,'rt') as rh:
        raxml_heads = rh.readlines()
        raxml_heads = [ll.strip() for ll in raxml_heads]
    with open(args.input_ncbiH,'rt') as nh:
        ncbi_heads = nh.readlines()
        ncbi_heads = [ll.strip() for ll in ncbi_heads]

    input_file = args.input_tree
    output_file = args.out_tree
    
    with open(output_file,'wt') as out_hdl:

        with open(input_file,'rt') as input_hdl:
            tree = input_hdl.read()
            for i,raxmlH in enumerate(raxml_heads):
                tree = tree.replace(raxmlH, ncbi_heads[i].split(" ")[0].replace(".","_"))

        out_hdl.write(tree)

    print("newick tree file with updated tip labels successfully saved")


if __name__ == "__main__":
    main()
