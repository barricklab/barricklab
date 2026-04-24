#!/usr/bin/env python3
"""
Subset reads from a fastq file based on length
"""

__author__ = "Ira Zibbu"
__version__ = "0.1.0"

""" imports """
import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Subset reads from a fastq file based on length')
parser.add_argument('--file', required=True, help='Path to fastq file')
parser.add_argument('--max', type=int, required=True, help='Maximum length of read to be subset')
parser.add_argument('--min', type=int, required=True, help='Minimum length of read to be subset')
parser.add_argument('--output', required=True, help='Name of file to output subset of reads')


def main(file, max_len, min_len, output):
    if not os.path.exists(file):
        raise IOError(f"Input file '{file}' does not exist.")

    filtered_reads = []

    for read in SeqIO.parse(file, "fastq"):
        if min_len <= len(read.seq) <= max_len:
            filtered_reads.append(read)

    SeqIO.write(filtered_reads, output, "fastq")


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.file, args.max, args.min, args.output)
