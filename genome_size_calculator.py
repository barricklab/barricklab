#!/usr/bin/env python3
"""
This script accepts a standard .gbk file, .fasta file or .gff3 file

"""

__author__ = "Ira Zibbu"
__version__ = "0.1.0"

""" imports """

from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
import logging
from collections import deque


""" Accept arguments """

parser = argparse.ArgumentParser(description='Command line code to compute length of the gbk, fasta or gff file. Right now only works for those with a single contig/ entry.')
parser.add_argument('--file', help='Path to file. Must end with .fasta, .fna, .gbk or .gff3')


def main(file):

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger('genome_size_calculator')
    length = 0
    if file.endswith("gbk"):
        try:
            record = SeqIO.read(file, "gb")
            length = len(record.seq)
        except Exception as e:
            logger.error(e)

    elif file.endswith("fasta") or file.endswith("fna"):
        try:
            record = SeqIO.read(file, "fasta")
            length = len(record.seq)
        except Exception as e:
            logger.error(e)

    elif file.endswith("gff3"):
        record_deque = deque() # store lines of the sequences as elements in a deque. Appending is O(1) in a deque
        is_sequence=False # will be false until code encounters the sequence data, which usually begins with a ">" fasta style header

        with open(file,'r') as sequence:
            for line in sequence:
                line = line.strip()
                if line[0] == ">":
                    is_sequence = True
                    continue
                if is_sequence:
                    record_deque.append(line)
            record = ''.join(record_deque)
            length = len(record)

    else:
        raise ValueError(f'Format of supplied file {file} could not be determined. Please ensure that the file ends in a valid extension')

    print(f"Length of file: {length} bp")

if __name__ == "__main__":
	args = parser.parse_args()
	main(args.file)
