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
"consensus_reads_with_stacks"
    no longer making distinction of SSCS/DCS
    using ustacks V1.48
    removed lambda functions from consensus base call function (replaced with simpler .itervalues and better use of key argument in max)
    reports updated
        gone SSCS/DCS/singletons
        added discrepancy counter (ie how many times was there disagreement among raw reads at a single base)
              non-consensus (ie how often could a consensus base not be identified
@author: ded
"""

import re
import argparse
import sys
import os
import subprocess

parser = argparse.ArgumentParser(description="Generate consensus sequences based on molecular indexes")
parser.add_argument("-c", "--constant_region", default="CAGT", help="expected constant region sequence")
parser.add_argument("-l", "--length", default=int(12), help="length of random MI region")
parser.add_argument("-f1", "--fastq1", help="fastq read1 file to check", required=True)
parser.add_argument("-f2", "--fastq2", help="fastq read2 file to check", required=True)
parser.add_argument("-cf", "--cutoff_frequency", default=float(0.5), type=float, help="minimum cutoff frequency required to call consensus base")
parser.add_argument("-p", "--prefix", help="prefix for output files", required=True)
parser.add_argument("-m", "--minimum_reads", default=2, type=int, help="minimum number of raw reads needed generate a consensus sequence")
parser.add_argument("--testing", action='store_true', default=False, help="break fastq read in after certain number of reads to increase speed")
# parser.add_argument("-v", "--verbose", action='store_true', help="increase output (not set to do anything)")  # TODO implement log function and verbose printing options within

# show help if no arguments passed
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Initialize key dictionary and search string, Verify proper command line arguments
main_dict = {"reads": [], "raw reads": 0, "dual MI reads": 0, "consensus sequence reads": 0, "consensus sequences": 0, "discrepancies removed": 0, "unknown consensus base": 0, "consensus": str()}  # Reads, counters, and consensus name
molecular_index_search_string = "^[ACTGN]{%d}%s" % (int(args.length), args.constant_region)
assert re.match("^[A-Za-z0-9_]+$", args.prefix), "\nNon-alphanumeric characters found:\n%s\nPlease restrict to letters, numbers, or underscore characters in prefix name" % [x for x in args.prefix if not re.match("[A-Za-z0-9_]", x)]
assert 0.5 <= args.cutoff_frequency <= 1.0, "Requiring a sub 50% frequency for consensus can lead to equal representation of 2 different bases. Requiring a frequency above 100% would mean no consenus could be generated. Please select a cutoff_frequency between 0.5 and 1"


def consensus_seq_gen(raw_sequences, consesnsus_freq=args.cutoff_frequency):
    """function to determine consensus sequence with given variable values. take in list of sequences, and consensus sequence above threshold."""
    assert len(list(set([len(x) for x in raw_sequences]))) == 1, "Unequal length sequences passed to consensus sequence generator.\n%s\n%s" % ([len(x) for x in raw_sequences], raw_sequences)
    consensus_read = []

    for i in xrange(len(raw_sequences[0])):  # i is the nucleotide index within the read
        consensus_dict = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0}  # set dictionary to empty for each base within the MI group
        for j in xrange(len(raw_sequences)):  # j is the read index within raw_sequences
            consensus_dict[raw_sequences[j][i]] += 1.0
        # if consensus_dict[max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key]))] / float(len(raw_sequences)) > consesnsus_freq:
        if max(consensus_dict.itervalues()) / float(len(raw_sequences)) > consesnsus_freq:  # is (frequency of highest represented base) > (frequency required)?
            # consensus_read.append(max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key])))
            consensus_read.append(max(consensus_dict, key=consensus_dict.get))  # append highest represented base to consensus read.
            main_dict["discrepancies removed"] += len(raw_sequences) - max(consensus_dict.itervalues())  # add number of reads with disagreeing base calls
        else:
            consensus_read.append("N")  # Most frequently observed base doesn't exceed consensus threshold
            main_dict["unknown consensus base"] += 1  # "True" base couldn't be determined therefore, count how often this occurs
    consensus_read = "".join(consensus_read)
    return consensus_read


def consensus_assignment():
    """generate consesnsus sequences and print them based on criteria"""
    if len(main_dict["reads"]) >= args.minimum_reads:
        consensus = consensus_seq_gen(main_dict["reads"])
        main_dict["consensus sequences"] += 1
        main_dict["consensus sequence reads"] += len(main_dict["reads"])
        print>>consensus_out, "\n".join(map(str, ["@_" + args.prefix.replace("_", "and") + "_consensus_" + main_dict["consensus"] + "_" + str(len(main_dict["reads"])) + "_" + str(main_dict["consensus sequences"]), consensus, "+", "I" * len(consensus)]))  # @_prefix_"consensus"_consensus.seq_#reads_#consensus.seq


def tag_parser():
    """function designed to parse .tags.tsv output from ustacks and print consensus sequences to *.concensus.fastq"""
    with open(args.prefix + ".tags.tsv", "r") as tags_in:
        for line in tags_in:
            line = line.rstrip().split("\t")
            if re.match("model", line[6]) or re.match("^#", line[0]):
                continue  # not doing anything with model lines, and header line at top starts with "#" in ustacks V 1.48
            elif re.match("consensus", line[6]):
                consensus_assignment()  # consesnsus sequence is first thing encountered, first time will correctly do nothing, but must be called after completion of the loop

                main_dict["reads"] = []  # reset sequences to empty lists
                main_dict["consensus"] = line[9]  # name will be used in header line of consensus.fastq file

            elif re.match("primary", line[6]):
                main_dict["reads"].append(line[8].split("_")[-1])  # non-MI sequence is after _ in fasta name line

            elif re.match("secondary", line[6]):
                assert False, "Secondary read alignment found, ustacks, not behaving as intended.\n%s" % line
            else:
                assert False, "Unknown sequence type identified, don't know how to handle this.\n%s\n%s" % (line[6], line)
        consensus_assignment()  # must run on last consensus sequence


def paired_end_fastq_read_in(fastq1, fastq2):
    """Function designed to read in paired end fastq files simultaneously, and pass each read pair to function to check for molecular index"""
    line_count = 0
    with open(fastq1, "r") as Fastq1, open(fastq2, "r") as Fastq2: #, open(args.prefix + ".log.txt", "w", 0) as log:
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
                # not doing anything with quality scores
                prep_for_stacks(header1, read1, header2, read2)
            if line_count % 200000 == 0:
                print>>log, line_count / 4, "reads processed"
                if args.testing:
                    break  # for testing subset of reads rather than full read list
    main_dict["raw reads"] = line_count / 4


def prep_for_stacks(header_1, read_1, header_2, read_2):
    """Steps: 1 Take in paired end reads, 2 Identify correct MI placement, 3 Output to new fasta file"""
    if re.match(molecular_index_search_string, read_1) and re.match(molecular_index_search_string, read_2):  # only work with read pairs that have correct MI on both reads
        # TODO slices are hardcoded. should be based on argparse variables
        f_tag = ">F-" + "-".join(header_1.split(":")[3:7]) + "_" + read_1[17:]
        f_read = read_1[:12] + read_2[:12]
        r_tag = ">R-" + "-".join(header_2.split(":")[3:7]) + "_" + read_2[17:]
        r_read = read_2[:12] + read_1[:12]
        print>>fasta_out, "\n".join(map(str, [f_tag, f_read, r_tag, r_read]))
        main_dict["dual MI reads"] += 1

# Read Fastq files in, output fasta files
assert not os.path.isfile(args.prefix + ".fasta"), "fasta file already exists. Please remove or rename existing file: %s.fasta" % args.prefix
with open(args.prefix + ".fasta", "w") as fasta_out, open(args.prefix + ".log.txt", "w", 0) as log:
    paired_end_fastq_read_in(args.fastq1, args.fastq2)

# Group reads by molecular index
stacks_command = ["ustacks", "-t", "fasta", "-f", str(args.prefix) + ".fasta", "-m", "1", "-M", "1", "-N", "0", "-H", "-p", "16", "--keep_high_cov", "-i", "1"]
external_stacks_call = subprocess.Popen(stacks_command, stdout=open(args.prefix + ".log.txt", "a"), stderr=open(args.prefix + ".log.txt", "a"))
external_stacks_call.wait()

# Read in molecular inedex, output consensus fastq file
assert not os.path.isfile(args.prefix + ".consensus.fastq"), "consensus file already exists. please remove or rename existing file: %s.consensus.fastq" % args.prefix
with open(args.prefix + ".consensus.fastq", "w") as consensus_out:
    tag_parser()

# Log final stats
with open(args.prefix + ".log.txt", "a") as log:
    print>>log, "\t".join(map(str, ["Sample", "raw reads", "dual MI reads", "consensus sequence reads", "consensus sequences", "discrepancies removed", "unknown consensus base"]))
    to_print = [args.prefix]
    for entry in ["raw reads", "dual MI reads", "consensus sequence reads", "consensus sequences", "discrepancies removed", "unknown consensus base"]:
        to_print.append(main_dict[entry])
    print>>log, "\t".join(map(str, to_print))
