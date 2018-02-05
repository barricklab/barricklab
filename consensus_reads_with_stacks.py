#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 11:15:09 2014
Stacks duplex seq

steps:
    1. read in R1 R2 fastq files together
    2. remove reads with incorrect MI
    3. output fasta file with following format: #unique number to come from header line, and be the coordinates of the read
        >F-readposition-R1(MI trimmed off)
        R1MIR2MI
        >R-unique#-R2(MI trimmed off)
        R2MIR1MI
    4. run stacks with following considerations:
        1. allow 1 MI to start a stack // -m 1
        2. allow 1 mismatch between stacks (corresponds to 3-4% error rate depending on if you count constant region 4bp) // -M 1
        3. first consideration eliminates the possible existence of secondary reads, therefore disable their mapping // -N 0
        4. happlotypes do not exist within MI, therefore disable // -H
        5. this assumed to be run on stampede, use maximal processors. currently set to 16 potential sight of softcode/dynamic choice // -p 16
        6. Keep highly represented sequences // --keep_high_cov
        7. Each sample required to have unique ID, as only running single sample at time, force all to standard 1 // -i 1  #TODO consider larger script to deal with multiples
    5. read in .tags.tsv file
        1. primary sequence names represent F/R reads
        2. store in  dictionary of dict ={F:[read name, read name, etc], R: [read name, read name]}
        3. dictionary resets list of sequences on consensus line,
        4. on each consensus, if length of either F or R list is greater than 1, generate consensus sequence.
    6. write output statistics to global log file
        1. figure out how to access singe file from multiple processes. need check to make sure is quickly appending and then closing the file and maybe some kind of wait/repeat if something returns error?
for gdname in 01_Data/*.gd;do name=$(echo $gdname|sed 's/.gd//'|sed 's/01_Data\///');echo "consensus_reads_with_stacks.py -f1 02_Downloads/$name*_R1_* -f2 02_Downloads/$name*_R2_* -p $name";done > commands
version:
2->3
    changes to consensus reads, not sscs/dcs. This means that all reads with same MI (regardless of R1/R2) will be combined to generate consensus
    need determination of when this should kick in (ie at the begining before stacks, or
@author: ded
"""

import re
import argparse
import sys
import os
import subprocess

parser = argparse.ArgumentParser(description="calculate how many reads have expected duplex seq tag")
parser.add_argument("-c", "--constant_region", default="CAGT", help="expected constant region sequence")
parser.add_argument("-l", "--length", default=int(12), help="length of random MI region")
parser.add_argument("-f1", "--fastq1", help="fastq read1 file to check")
parser.add_argument("-f2", "--fastq2", help="fastq read2 file to check")
parser.add_argument("-cf", "--cutoff_frequency", default=float(0.5), type=float, help="minimum cutoff frequency required to call consensus base")
parser.add_argument("-p", "--prefix", help="prefix for output files")
parser.add_argument("-m", "--minimum_reads", default=2, type=int, help="minimum number of reads needed to support SSCS reads")
# parser.add_argument("-v", "--verbose", action='store_true', help="increase output (not set to do anything)")
# parser.add_argument("-i", "--intermediate_fastqs", action='store_true', help="remove MI and output 'double, single, zero'.fastq files based on presence of 2,1,0 MI. DOES NOT EFFECT DCS, SSCS, or Singleton output files.")  # depreciated
parser.add_argument("-d", "--DCS", action='store_true', default=False, help="calculate DCS sequence, off by default. SSCS on as well.")  # TODO set SSCS to be separate from DCS, and for SSCS to be on by default
parser.add_argument("--testing", action='store_true', default=False, help="break fastq read in after certain number of reads to increase speed")
# parser.add_argument("-s", "--SSCS", action='store_true', default = False, help="calculate SSCS sequence, off by default. IF DCS specified, automatically on")  # TODO reimplement for general ability to trim, sscs, or dcs more fully.
# parser.add_argument("-t", "--trimmed", help="enter 2 for requiring correctly identified MI on both reads, enter 1 for requiring correct MI on single read. non")#hasto be on for anything to happen?

#show help if no arguments passed
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()


main_dict = {"F": [], "R": [], "reads": [], "total sequences": 0, "total consensus sequences": 0, "consensus sequences": 0, "DCS sequences": 0, "SSCS sequences": 0, "double singleton sequences": 0, "consensus": str()}  # F, R sequences, counters, and consensus name, double singleton sequences counter to be used for determining value of 2x read sscs
molecular_index_search_string = "^[ACTGN]{%d}%s" % (int(args.length), args.constant_region)
assert re.match("^[A-Za-z0-9_]+$", args.prefix), "\nNon-alphanumeric characters found:\n%s\nPlease restrict to letters, numbers, or underscore characters in prefix name" % [x for x in args.prefix if not re.match("[A-Za-z0-9_]", x)]


def consensus_seq_gen(raw_sequences, consesnsus_freq=args.cutoff_frequency):
    """function to determine consensus sequence with given variable values"""
    """take in list of sequences, and consensus sequence above threshold.
    Largely copied from Scott Kennedy's consensusMaker.py script used by Schmitt et all
    partially copied from Scott Kennedy's script
    ok, not really copied from tht script anymore"""
    assert len(list(set([len(x) for x in raw_sequences]))) == 1, "Unequal length sequences passed to consensus sequence generator.\n%s\n%s" % ([len(x) for x in raw_sequences], raw_sequences)
    consensus_read = []

    for i in xrange(len(raw_sequences[0])):  # Count the types of nucleotides at a position in a read. i is the nucleotide index within a read in list_of_sequences
        consensus_dict = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0}  # set dictionary  to empty for each base within the readgroup
        for j in xrange(len(raw_sequences)):  # Do this for every read that comprises a SMI group. j is the read index within list_of_sequences
            consensus_dict[raw_sequences[j][i]] += 1.0
        if consensus_dict[max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key]))] / float(len(raw_sequences)) > consesnsus_freq:  # this would fall apart if a cutoff of less than 0.5 is used
            consensus_read.append(max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key])))
        else:
            consensus_read.append("N")  # previous version of script had a counter running for how often consensus could not be reached
    consensus_read = "".join(consensus_read)
    return consensus_read


def consensus_assignment():
    """generate consesnsus sequences and print them based on criteria"""
    # with open(args.prefix + ".sscs.fastq", "a") as SSCS_out:
    # F = main_dict["F"] ## depreciated on v3
    # R = main_dict["R"] ## depreciated on v3
    if len(main_dict["reads"]) >= args.minimum_reads:
        consensus = consensus_seq_gen(main_dict["reads"])
        main_dict["consensus sequences"] += 1
        print>>consensus_out, "\n".join(map(str, ["@_" + args.prefix.replace("_", "and") + "_consensus_" + main_dict["consensus"] + "_" + str(len(main_dict["reads"])) + "_" + str(main_dict["consensus sequences"]), consensus, "+", "I" * len(consensus)]))  # @_prefix_"consensus"_consensus.seq_#reads_#consensus.seq
    # if len(F) >= args.minimum_reads and len(R) >= args.minimum_reads:  # both F and R can generate SSCS
    #     main_dict["DCS sequences"] += 1
    #     main_dict["SSCS sequences"] += 2
    #     F_consensus = consensus_seq_gen(F)
    #     R_consensus = consensus_seq_gen(R)
    #     if args.DCS:
    #         # with open(args.prefix + ".dcs.fastq", "a") as DCS_out:
    #         DCS_consensus = consensus_seq_gen(F_consensus + R_consensus, 0.5)
    #         print>>DCS_out, "\n".join(map(str, ["@_" + args.prefix + "_DCS_" + main_dict["consensus"] + "_" + str(main_dict["DCS sequences"]), DCS_consensus, "+", "I" * len(DCS_consensus)]))
    #     print>>SSCS_out, "\n".join(map(str, ["@_" + args.prefix + "_SSCS_" + main_dict["consensus"] + "_" + str(main_dict["SSCS sequences"] - 1), F_consensus, "+", "5" * len(F_consensus)]))  # -1 on the counter to account for this being a DCS qualifying read, and wanting SSCS numbering to be separate
    #     print>>SSCS_out, "\n".join(map(str, ["@_" + args.prefix + "_SSCS_" + main_dict["consensus"] + "_" + str(main_dict["SSCS sequences"]), R_consensus, "+", "5" * len(R_consensus)]))
    # elif len(F) >= args.minimum_reads:  # only F can generate SSCS
    #     main_dict["SSCS sequences"] += 1
    #     F_consensus = consensus_seq_gen(F)
    #     print>>SSCS_out, "\n".join(map(str, ["@_" + args.prefix + "_SSCS_" + main_dict["consensus"] + "_" + str(main_dict["SSCS sequences"]), F_consensus, "+", "5" * len(F_consensus)]))
    # elif len(R) >= args.minimum_reads:  # only R can generate SSCS
    #     main_dict["SSCS sequences"] += 1
    #     R_consensus = consensus_seq_gen(R)
    #     print>>SSCS_out, "\n".join(map(str, ["@_" + args.prefix + "_SSCS_" + main_dict["consensus"] + "_" + str(main_dict["SSCS sequences"]), R_consensus, "+", "5" * len(R_consensus)]))
    # else:  # neither F or R can generate SSCS
    #     if len(F) == 1 and len(R) == 1:
    #         main_dict["double singleton sequences"] += 1
    #         pass  # TODO implement singletons here?


def tag_parser():
    """function designed to parse .tags.tsv output from ustacks and print consensus sequences to fil eas either DCS or SSCS"""
    # assert not os.path.isfile(args.prefix + ".sscs.fastq"), "SSCS file already exists. please move remove or rename existing file: %s.sscs.fastq" % args.prefix
    with open(args.prefix + ".tags.tsv", "r") as tags_in:
        # if args.DCS:
            # assert not os.path.isfile(args.prefix + ".dcs.fastq"), "DCS selected, but DCS file already exists, please move, remove, or rename existing file %s.dcs.fastq" % args.prefix
        for line in tags_in:
            line = line.rstrip().split("\t")
            if re.match("^#", line[0]):
                continue  # ustacks V 1.48 has header line at top
            if re.match("model", line[6]):
                continue  # not doing anything with model lines
            elif re.match("consensus", line[6]):
                consensus_assignment()  # consesnsus sequence is first thing encountered, first time through will do nothing, this is expected, and hence the call after completion of the loop
                # main_dict["F"] = []  # reset sequences to empty lists ##depreciated on v3
                # main_dict["R"] = []  # reset sequences to empty lists ##depreciated on v3
                main_dict["reads"] = []  # reset sequences to empty lists
                main_dict["consensus"] = line[9]
                main_dict["total consensus sequences"] += 1
            elif re.match("primary", line[6]):
                # main_dict[line[8].split("-")[0].replace(">", "")].append(line[8].split("_")[-1]) ##depreciated on v3
                main_dict["total sequences"] += 1
                main_dict["reads"].append(line[8].split("_")[-1])
                if args.testing and main_dict["total sequences"] % 200000 == 0:
                    break  # for testing subset rather than full read list
            elif re.match("secondary", line[6]):
                assert False, "Secondary read alignment found, ustacks, not behaving as intended.\n%s" % line
            else:
                assert False, "Unknown sequence type identified, don't know how to handle this.\n%s\n%s" % (line[6], line)
        consensus_assignment()  # must run on last consensus sequence


def paired_end_fastq_read_in(fastq1, fastq2):
    """Function designed to read in paired end fastq files simultaneously, and do something with them
    Has the additional goal of making this something that can be loaded/copied into other scripts"""
    line_count = 0
    with open(fastq1, "r") as Fastq1, open(fastq2, "r") as Fastq2, open(args.prefix + ".log.txt", "w", 0) as log:
        for line1 in Fastq1:
            line1 = line1.rstrip()
            line2 = Fastq2.next().rstrip()
            line_count += 1
            if line_count % 4 == 1:  # line = header
                header1 = line1
                header2 = line2
                assert (header1.split(" ")[0] == header2.split(" ")[0]) and (len(header1.split(" ")) == len(header1.split(" "))) and len(header1.split(" ")) == 2, "Header lines of R1 and R2 do not agree with expectations: single space in header line, and header before space is identical for both reads.\nHeaders:\n%s\n%s" % (header1, header2)  # this may be inaccurate for different Illumina outputs
            elif line_count % 4 == 2:  # line = read
                read1 = line1
                read2 = line2
                assert re.match('^[ACTGN]*$', read1) and re.match('^[ACTGN]*$', read2), "Non-sequence characters found on one of th following sequence lines:\n%s\n%s" % (line1, line2)
            elif line_count % 4 == 3:  # line = optional sequence descriptor. assume line to use 'standard' "+" notation rather than actual descriptor.
                assert line1 == line2 == "+", "line1 or line2 does not display expected '+' sign on line3 here:\nline1: %s\nline2: %s" % (line1, line2)
            elif line_count % 4 == 0:  # line = quality
                quality1 = line1
                quality2 = line2
                # TODO this is where the function should actually do something with the read pair, not sure how to make this generizable to other functions ... perhaps chain function calls so the paired_end_fastq_read_in function takes a 3rd parameter of another function ... should pass header, read, and quality out, and only use appropriate
                prep_for_stacks(header1, read1, header2, read2)
            if line_count % 200000 == 0:
                print>>log, line_count / 4, "reads processed"
                if args.testing:
                    break  # for testing subset of reads rather than full read list


def prep_for_stacks(header_1, read_1, header_2, read_2):
    """Steps: 1 Take in paired end reads, 2 Identify correct MI placement, 3 Output to new fasta file"""
    if re.match(molecular_index_search_string, read_1) and re.match(molecular_index_search_string, read_2):  # only work with read pairs that have correct MI on both reads
        f_tag = ">F-" + "-".join(header_1.split(":")[3:7]) + "_" + read_1[17:]
        f_read = read_1[:12] + read_2[:12]
        r_tag = ">R-" + "-".join(header_2.split(":")[3:7]) + "_" + read_2[17:]
        r_read = read_2[:12] + read_1[:12]
        print>>fasta_out, "\n".join(map(str, [f_tag, f_read, r_tag, r_read]))

assert not os.path.isfile(args.prefix + ".fasta"), "fasta file already found"  # TODO set this up as intermediate starting point?
with open(args.prefix + ".fasta", "w") as fasta_out:
    paired_end_fastq_read_in(args.fastq1, args.fastq2)

stacks_command = ["ustacks", "-t", "fasta", "-f", str(args.prefix) + ".fasta", "-m", "1", "-M", "1", "-N", "0", "-H", "-p", "16", "--keep_high_cov", "-i", "1"]
external_stacks_call = subprocess.Popen(stacks_command, stdout=open(args.prefix + ".log.txt", "a"), stderr=open(args.prefix + ".log.txt", "a"))
external_stacks_call.wait()

assert not os.path.isfile(args.prefix + ".consensus.fastq"), "consensus file already exists. please remove or rename existing file: %s.consensus.fastq" % args.prefix
with open(args.prefix + ".consensus.fastq", "w") as consensus_out:
    tag_parser()

# assert not os.path.isfile(args.prefix + ".sscs.fastq"), "SSCS file already exists. please move remove or rename existing file: %s.sscs.fastq" % args.prefix
# if args.DCS:
#     assert not os.path.isfile(args.prefix + ".dcs.fastq"), "DCS selected, but DCS file already exists, please move, remove, or rename existing file %s.dcs.fastq" % args.prefix
#     with open(args.prefix + ".dcs.fastq", "w") as DCS_out, open(args.prefix + ".sscs.fastq", "a") as SSCS_out:
#         tag_parser()
# else:
#     with open(args.prefix + ".sscs.fastq", "w") as SSCS_out:
#         tag_parser()

with open(args.prefix + ".log.txt", "a") as log:
    for entry in ["total sequences", "total consensus sequences", "consensus sequences", "DCS sequences", "SSCS sequences", "double singleton sequences"]:
        print>>log, entry, "\t", main_dict[entry]
