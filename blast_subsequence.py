#!/usr/bin/env python3
"""
This script accepts a reference FASTA file, a range <int>-<int>, a query FASTA file. It extracts a subsequence from the query, and uses that to blast against the query FASTA.
Requires a conda env with blast and biopython installed

"""

__author__ = "Ira Zibbu"
__version__ = "0.1.0"

""" imports """

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from Bio.SeqRecord import SeqRecord
import logging
import os
import subprocess


""" Accept arguments """

parser = argparse.ArgumentParser(description='This script accepts a reference FASTA file, a range <int>-<int> (inclusive to start, excludes stop) a query FASTA file. It extracts a subsequence from the query, and uses that to blast against the query FASTA. Requires a conda env with blast installed')
parser.add_argument('--reference', help='Path to reference FASTA file')
parser.add_argument('--start', help='<int> 1-based index of the start of the range')
parser.add_argument('--stop', help='<int> 1-based index of the stop of the range')
parser.add_argument('--query', help='Path to query FASTA file')
parser.add_argument('--output', help='Path to txt file')

def read_fasta_files(reference_path,query_path):

    """
    Parse the reference and query fasta file with biopython

    Args:
        reference_path: path to the reference fasta file
        query_path: path to the query fasta file

    Returns:
        reference: biopython seqrecord for the reference file
        query: biopython seqrecord for the query file

    Raises:
        IOError: reference/query path is incorrect or file does not exist at specified path
        Generic exception:

    """

    if not(os.path.exists(reference_path)):
        raise IOError(f"{reference_path} reference file does not exist")

    if not(os.path.exists(query_path)):
        raise IOError(f"{query_path} reference file does not exist")

    reference = None
    query = None

    try:
        reference = SeqIO.read(reference_path, "fasta")
    except Exception as e:
        print(f"Error parsing reference file: {reference_path}. Exception: {e}")
        raise

    try:
        query = SeqIO.read(query_path, "fasta")
    except Exception as e:
        print(f"Error parsing query file: {query_path}. Exception: {e}")
        raise

    if reference is None or query is None:
        raise ValueError("Both reference and query files must contain valid FASTA records.")

    logging.info('Reference and query file read successfully')
    return reference, query

def extract_subsequence(reference,start,stop):

    """
    Extract the subsequence start:stop from the reference sequence. Writes out tmp.fasta in the wd which corresponds to this subsequence

    Args:
        reference: seqrecord representing the referenece
        start: 1-based index start
        stop: 1-based index stop (excluding)

    """

    # shift to 1-based indexing
    start=int(start)+1
    stop=int(stop)+1

    subsequence = SeqRecord(
        reference.seq[start:stop],
        id=reference.id,
        description=f"{reference.id} bases {start} to {stop}",
    )

    with open("tmp.fasta", "w") as output_handle:
        SeqIO.write(subsequence, output_handle, "fasta")
        logging.info('tmp.fasta created')

    return subsequence


def run_blast(query_path,output):

    """
    Use blastn to align the subsequence to the reference

    Args:
        reference_path = path to reference fasta file
        output = file name to store blast output

    Returns:
        None

    Raises:
        CalledProcessError: An error while running the blast subprocess

    """

    blast_command = [
        'blastn','-subject',query_path,'-query','tmp.fasta']

    try:
        # Run the subprocess
        with open(output, 'w') as output_file:
            output_file.write(f"Query of blast is: {query_path}")
            output_file.write(f"Subject of blast is: {subject_path}")
            result = subprocess.run(blast_command, stdout=output_file)
            logging.info("Subprocess blast completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess blast: {e}")

def test():

    reference_path="test/Ara+2_10000gen_4531A.fasta"
    query_path="test/REL606.fasta"

    reference,query = read_fasta_files(reference_path,query_path)
    assert reference.id == "REL606"

    subsequence = extract_subsequence(reference,0,5)

    assert subsequence.seq == "GCTTT"


def main(reference_path,query_path,start,stop,output):

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger('blast_subsequence')

    #test()
    reference,query = read_fasta_files(reference_path,query_path)
    extract_subsequence(reference,start,stop)
    run_blast(query_path,output)

    os.remove("tmp.fasta")
    logging.info("Cleaning up tmp files")


if __name__ == "__main__":
	args = parser.parse_args()
	main(args.reference,args.query,args.start,args.stop,args.output)
