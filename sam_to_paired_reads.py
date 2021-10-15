#!/usr/bin/env python

# -*- coding: utf-8 -*-

import os
import re
import argparse
import sys
import re
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser(description="Extract all read pairs appearing in SAM file")
parser.add_argument("-s", "--sam", help="SAM file path.", required=True)
parser.add_argument("-I", "--fastq1", help="Input FASTQ path for 1st read in each pair", required=True)
parser.add_argument("-i", "--fastq2", help="Input FASTQ path for 2nd read in each pair", required=True)
parser.add_argument("-O", "--output1", help="Output FASTQ path for 1st read in each pair", required=True)
parser.add_argument("-o", "--output2", help="Output FASTQ path for 1st read in each pair", required=True)

# show help if no arguments passed
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

print(args)

samFile = open(args.sam, "r")

fastq1File = open(args.fastq1, "r")
fastq2File = open(args.fastq2, "r")

output1File = open(args.output1, "w")
output2File = open(args.output2, "w")

seqIDs = {}
for line in samFile:
    seqIDs[line.split("\t")[0]] = True
    
print(seqIDs)

fastq1it = FastqGeneralIterator(fastq1File)
fastq2it = FastqGeneralIterator(fastq2File)
for (id1, seq1, qual1) in fastq1it:
    #print(id1)
    (id2, seq2, qual2) = next(fastq2it)
    #print(id2)
    if id1 in seqIDs or id2 in seqIDs:
        output1File.write("@%s\n%s\n+\n%s\n" % (id1, seq1, qual1))
        output2File.write("@%s\n%s\n+\n%s\n" % (id2, seq2, qual2))
