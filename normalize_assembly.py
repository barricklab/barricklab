#! /usr/bin/env python3

# This script:
#   Renames contigs to match or be in order 
#   Sorts contigs in reverse order by length
#   Reindexes (shifts) the sequence in the input file to a requested subsequence

# ---- Imports ----
import argparse
import os
import importlib.util
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def reindex_sequence(in_sequence, in_sequence_name, in_reindex_bases):

    print("reindex")
    out_sequence = in_sequence

    found = False
    if (not found): 
        search_sequence = in_sequence + in_sequence[0:len(in_reindex_bases)]
        try: # the forward sequence
            reindex_start = search_sequence.index(in_reindex_bases)
            out_sequence = in_sequence[reindex_start:]+in_sequence[:reindex_start]
            print("Reindexed {} contig to sequence {}".format(in_sequence_name, in_reindex_bases))
            found = True
        except:
            pass
        
    if (not found): 
        contig_rc = Seq(in_sequence).reverse_complement()
        search_sequence = contig_rc + contig_rc[0:len(in_reindex_bases)]
        try:
            reindex_start = search_sequence.index(in_reindex_bases)
            out_sequence = contig_rc[reindex_start:]+contig_rc[:reindex_start]
            print("Reindexed {} contig to sequence {} as reverse complement".format(in_sequence_name, in_reindex_bases))
            found = True
        except:
            pass


    # Abandon hope if we can't find the reindex bases
    if (not found):
        print("Could not reindex {} to bases {}".format(in_sequence_name, in_reindex_bases))

    return(out_sequence)

def main():
    # Get arguments
    parser = argparse.ArgumentParser(description="Reindex sequence file to a conserved sequence")
    parser.add_argument("-r", "--reference", help="Bases to reindex on", action="store", default=None)
    parser.add_argument("-i", "--input", help="Folder containing target files", action="store", required=True)
    parser.add_argument("-o", "--output", help="Folder containing target files", action="store", required=True)
    parser.add_argument("-t", "--file-type", help="Filetype to reindex [fasta]", action="store", default="fasta")
    parser.add_argument("-s", "--sort", help="Sort contigs, longest to shortest. Sorting happens before renaming", action='store_true')
    parser.add_argument("-c", "--copy-names", help="Rename contig names to match the reference in order.", action='store_true')
    parser.add_argument("-n", "--new-names", help="Rename contigs with  this prefix followed by their index in the file", action="store")
    parser.add_argument("-x", "--reindex", help="Reindex (rotate) contigs so they start with the same bases as the reference", action='store_true')
    parser.add_argument("-y", "--reindex-initial-bases", help="Reindex (rotate) contigs in so they start with the same number of identical nucleotides to the reference", action="store", type=int, default=24)
    parser.add_argument("-z", "--reindex-sequence", help="Reindex (rotate) all contigs to begin with this sequence if possible", action="store")

    args = parser.parse_args()

    # Check if filetype is valid
    valid_file_types = ["fasta"]
    if args.file_type.lower() not in valid_file_types:
        print("Invalid file type. Valid file types are: {}".format(", ".join(valid_file_types)))
        exit(1)

    # Load reference file

    reference_seqs = []
    if (args.reference != None):
        print("Reference File: ", args.reference)
        for record in SeqIO.parse(args.reference, "fasta"):
            reference_seqs.append({'id' : record.id, 'seq' : record.seq})
        print("  * Contigs:" + str([d['id'] for d in reference_seqs]) )

    # Load input file
    input_seqs = []
    print("Input File: ", args.input)
    for record in SeqIO.parse(args.input, "fasta"):
        input_seqs.append({'id' : record.id, 'seq' : record.seq})
    print("  * Contigs:" + str([d['id'] for d in input_seqs]) )

    # Handle sorting
    if (args.sort):
        print("SORTING")
        input_seqs = sorted(input_seqs, key=lambda d: -len(d['seq'])) 

    # Handle renaming (copying)
    if (args.copy_names):
        print("COPYING NAMES FROM REFERENCE")
        assert args.reference != None, "Must provide reference for this operation."
        assert len(reference_seqs)==len(input_seqs), "Reference and input files must have the same number of sequences."
        for i, item in enumerate(input_seqs):
            item['id'] = reference_seqs[i]['id'] 
    
    # Handle renaming (by index)
    if (args.new_names):
        print("RENAMING BY INDEX")
        assert len(reference_seqs)==len(input_seqs), "Reference and input files must have the same number of sequences."
        for i, item in enumerate(input_seqs):
            item['id'] = args.new_names + str(i+1)

    # Reindex
    if (args.reindex):
        print("REINDEXING TO REFERENCE")
        assert args.reference != None, "Must provide reference for this operation."
        number_of_bases_to_match = args.reindex_initial_bases
        for i, item in enumerate(input_seqs):
            item['seq'] = reindex_sequence(item['seq'], item['id'], reference_seqs[i]['seq'][0:number_of_bases_to_match-1])

    if (args.reindex_sequence):
        print("REINDEXING TO SEQUENCE")
        for i, item in enumerate(input_seqs):
            item['seq'] = reindex_sequence(item['seq'], item['id'], args.reindex_sequence)

    #Write output file
    with open(args.output, "w") as output_handle:
        for i in input_seqs:
            simple_seq_r = SeqRecord(Seq(i['seq']))
            simple_seq_r.id = i['id']
            simple_seq_r.description = ""
            SeqIO.write(simple_seq_r, output_handle, "fasta")


if __name__ == '__main__':
    main()

