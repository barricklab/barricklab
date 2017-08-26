#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tuesday March 21 2017 5:07pm
program designed to read in plasmid sequence, and identify all sub-sequences repeated more than once.

work in progress
"""

import argparse
import re
from collections import defaultdict
import sys

parser = argparse.ArgumentParser(description="read in .fastq.gz file names, generate .gd files. Script allows for single '_' in sample names")
parser.add_argument("-i", "--input", help="Plasmid file in fasta format", required=True)
# parser.add_argument("-f", "--format", help="Format of plasmid file. Supported types: fasta, gbk, more", required=True)  ##TODO add other formats
parser.add_argument("-r", "--repeat_size", type=int, help="Minimum size of repeat to find. must be >=3", default="5")
parser.add_argument("-l", "--locations", action='store_true', help="Print all locations of repeats. Can create a large number of columns and make difficult/impossible to open in Excel.")
parser.add_argument("-o", "--output", default="new_gd_files", help="Output tsv file. Will default to 'plasmid_name + _repeats.tsv'")
# parser.add_argument("-s", "--sra", action='store_true', help="Files are on the NCBI SRA archive.")
# parser.add_argument("-v", "--verbose", action='store_true', help="Increase output")
# parser.add_argument("-m", "--meta", help="Optional meta data tsv file with headers. 'sample' required, any of the following allowed: 'population', 'time' (meaning generation), 'treatment', 'clone'")


# show help if no arguments passed
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# sequence read in
with open(args.input, "r") as f:
    header_line = True
    sequence = ''
    for line in f:
        line = line.rstrip().upper()
        if header_line:
            assert re.match(">", line), "File does not appear to be fasta format as first line does not begin with '>':\n%s" % line
            header_line = False
            continue
        assert re.match("^[ACTG]*$", line), "Non sequence character found:\n%s" % line
        sequence += line

repeat_dict = defaultdict(dict)

size = args.repeat_size

print len(sequence)
for i1, b1 in enumerate(sequence):
    match_found = True
    size = args.repeat_size
    number_loops = 0  # for testing if while loop is freezing
    while match_found:
        number_loops += 1
        # print number_loops  # for testing if while loop is freezing
        search_string = sequence[i1:i1 + size]
        if len(search_string) < size:  # at end of sequence start running out of room, can probably be eliminated by changing search space
            search_string += sequence[0:size - len(search_string)]
        assert len(search_string) == size

        ## block originally used to avoid searching for strings already searched, but generates additional key 1 larger than searched
        ## changed to only storing if new
        # if search_string in repeat_dict[size]:
        #     size += 1  # this was found multiple times before, therefore know it will match multiple times, therefore increase size
        #     continue  # have already searched for subsequence and found it multiple times, no need to research for it

        matches = [(m.start(), m.end()) for m in re.finditer(search_string, sequence + sequence[0: size - 1])]  # a list of 2 item tuple of start loc and end loc
        if len(matches) > 1:  # multiple matches found
            if search_string not in repeat_dict[size]:
                repeat_dict[size][search_string] = matches
            else:
                assert repeat_dict[size][search_string] == matches

            size += 1
            matches = []  # reset matches to empty list
        else:  # does not repeat
            match_found = False


# for size in repeat_dict:
#     for seq in repeat_dict[size]:
#         print seq, repeat_dict[size][seq]
#     #print size, len(repeat_dict[size])
# #print sequence
# for size in repeat_dict:
#     print size, "\t", len(repeat_dict[size])
# print sorted(repeat_dict.keys(), reverse=True)

no_overlaps = defaultdict(dict)

for size in sorted(repeat_dict.keys(), reverse=True):
    for r_seq in repeat_dict[size]:
        uniq_seq = True
        assert r_seq not in no_overlaps
        for n_o_seq in no_overlaps:
            if re.search(r_seq, n_o_seq):  # sequence is a subsequence of existing sequence
                if len(repeat_dict[size][r_seq]) == len(no_overlaps[n_o_seq]):  # if there are different numbers then interested in finding it for completeness
                    overlap_list = []
                    for repeat_loc in repeat_dict[size][r_seq]:  # loop through all locations of this sequence
                        assert repeat_loc[0] < repeat_loc[1], repeat_loc
                        overlapped_seq = False
                        for no_overlap_loc in no_overlaps[n_o_seq]:  # loop through all locations of the sequence that the test sequence is a subset of
                            if no_overlap_loc[0] <= repeat_loc[0] <= repeat_loc[1] <= no_overlap_loc[1]:
                                overlapped_seq = True
                                break  # repeat locations have overlapped no_overlap locations, move to next repeat location
                        overlap_list.append(overlapped_seq)
                    assert len(overlap_list) == len(repeat_dict[size][r_seq])
                    if len(set(overlap_list)) == 1 and overlap_list[0] == True:
                        uniq_seq = False
        if uniq_seq:
            no_overlaps[r_seq] = repeat_dict[size][r_seq]

final_dict = defaultdict(list)

for final_seq in no_overlaps:
    final_dict[len(final_seq)].append(final_seq)

print "\t".join(map(str, ["Size", "Number of Sequences", "Sequences"]))

for size in sorted(final_dict.keys(), reverse=True):
    to_print = [size, len(final_dict[size])]
    to_print += [x for x in final_dict[size]]
    print "\t".join(map(str, to_print))

headings = ["Size", "Sequence", "Number of Repeats"]
if args.locations:
    headings.append("Locations")
print "\t".join(map(str, headings))

for size in sorted(final_dict.keys(), reverse=True):
    for seq in final_dict[size]:
        toprint = [size, seq, len(no_overlaps[seq])]
        if args.locations:
            for loc in no_overlaps[seq]:
                x = "%s - %s" % (loc[0], loc[1])
                toprint.append(x)
        print "\t".join(map(str, toprint))

for seq in no_overlaps:
    assert len(re.findall(seq, sequence + sequence[0:len(seq)-1])) > 1, re.findall(seq, sequence + sequence[0:len(seq)-1])
    assert len(set(re.findall(seq, sequence + sequence[0:len(seq)-1]))) == 1, set(re.findall(seq, sequence + sequence[0:len(seq)-1]))
    for m in re.finditer(seq, sequence + sequence[0: len(seq) - 1]):
        assert m.start() < len(sequence)
