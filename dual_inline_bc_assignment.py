#!/opt/apps/intel15/python/2.7.9/bin/python
__author__ = 'ded'
"""
Script designed to read in dual internal barcode files, and output reads with internal barcodes stripped from reads.
"""

import re
import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Read paired end fastq files in, assign reads to different files based on barcodes")
parser.add_argument("-l", "--length", default=7, help="length of internal barcodes", type=int)
parser.add_argument("-f1", "--fastq1", help="fastq read1 file to check")
parser.add_argument("-f2", "--fastq2", help="fastq read2 file to check")
parser.add_argument("-e", "--expected", help="2 column tsv file, no headers of output file name, and expected barcode sequence. Barcode sequence should be R1R2 no spaces or punctuation.")
parser.add_argument("-c", "--combine", help="combine R1 and R2 into a single fastq output file", default=False, action='store_true')
parser.add_argument("-v", "--verbose", help="print progress of fastq read in", action='store_true')
args = parser.parse_args()

# Set up read dictionary
read_dict = defaultdict()
output_dict = {"unknown": "unknown_barcodes"}
read_stats = {"unknown": 0}
with open(args.expected, "r") as f:
    for line in f:
        line = line.rstrip().split("\t")
        assert len(line) == 2, "more than 2 columns detected %s" % line
        assert len(line[1]) == 2 * args.length, "barcode given (%s) is different than length (%i) of expected barcode" % (line[1], 2 * args.length)
        line[1] = "A" + line[1]  # this needs more validation
        assert line[1] not in read_dict, "identical barcodes given %s" % line[1]
        read_stats[line[1]] = 0
        if args.combine:
            read_dict[line[1]] = []
        else:
            read_dict[line[1]] = {"_R1": [], "_R2": []}
        output_dict[line[1]] = line[0].rstrip(".fastq")
if args.combine:
    read_dict["unknown"] = []
else:
    read_dict["unknown"] = {"_R1": [], "_R2": []}

# Verify output file names are valid
for entry in output_dict:
    if args.combine:
        assert not os.path.exsists(output_dict[entry] + ".fastq"), "%s file already exists. Move existing file to new location or delete" % (output_dict[entry] + ".fastq")
    else:
        for _ in ["_R1", "_R2"]:
            assert not os.path.exists(output_dict[entry] + _ + ".fastq"), "%s file already exists. Move existing file to new location or delete" % (output_dict[entry] + _ + ".fastq")

# Read fastq1 an fastq2 in together
line_count = 0
with open(args.fastq1, "r") as Fastq1, open(args.fastq2, "r") as Fastq2:
    for line1 in Fastq1:
        line1 = line1.rstrip()
        line2 = Fastq2.next().rstrip()  # advance read 2 file to keep in sync
        line_count += 1
        if line_count % 4 == 1:  # line = header
            header1 = line1
            header2 = line2
            assert (header1.split(" ")[0] == header2.split(" ")[0]) and (len(header1.split(" ")) == len(header1.split(" "))) and len(header1.split(" ")) == 2, "Header lines of R1 and R2 do not agree with expectations: single space in header line, and header before space is identical for both reads.\nHeaders:\n%s\n%s" % (header1, header2)  # this may be inaccurate for different illumina outputs
        elif line_count % 4 == 2:  # line = read
            assert re.match('^[ACTGN]*$', line1) and re.match('^[ACTGN]*$', line2), "Non-sequence characters found on one of th following sequence lines:\n%s\n%s" % (line1, line2)
            read1 = line1
            read2 = line2
        elif line_count % 4 == 3:  # line = optional sequence descriptor. assume line to use 'standard' "+" notation rather than actual descriptor.
            assert line1 == line2 == "+", "line1 or line2 does not display expected '+' sign on line3 here:\nline1: %s\nline2: %s" % (line1, line2)
        elif line_count % 4 == 0:  # line = quality score
            quality1 = line1[args.length + 1:]
            quality2 = line2[args.length:]
            if args.combine:
                try:
                    read_dict[read1[:args.length + 2] + read2[:args.length + 1]].extend([header1, read1[args.length + 1:], "+", quality1, header2, read2[args.length:], "+", quality2])
                except KeyError:  # barcodes not specified, therefore sore as unknown
                    #  NOTE that the typical "+" symbol on line 3 is changed to what the unidentified barcode was, with first half corresponding to read 1 and second half corresponding to read 2
                    read_dict["unknown"].extend([header1, read1[args.length + 1:], read1[:args.length + 2] + read2[:args.length + 1], quality1, header2, read2[args.length:], read1[:args.length + 2] + read2[:args.length + 1], quality2])
            else:  # R1 and R2 to be kept separate in final output
                try:
                    read_dict[read1[:args.length + 2] + read2[:args.length + 1]]["_R1"].extend([header1, read1[args.length + 1:], "+", quality1])
                    read_dict[read1[:args.length + 2] + read2[:args.length + 1]]["_R2"].extend([header2, read2[args.length:], "+", quality2])
                except KeyError:  # barcodes not specified, therefore sore as unknown
                    #  NOTE that the typical "+" symbol on line 3 is changed to what the unidentified barcode was, with first half corresponding to read 1 and second half corresponding to read 2
                    read_dict["unknown"]["_R1"].extend([header1, read1[args.length + 1:], read1[:args.length + 2] + read2[:args.length + 1], quality1])
                    read_dict["unknown"]["_R2"].extend([header2, read2[args.length:], read1[:args.length + 2] + read2[:args.length + 1], quality2])
        if args.verbose and line_count % 200000 == 0:
            print line_count / 4, "reads processed"
            # break  # Uncomment for testing subset of reads rather than full read list
        if line_count % 40000000 == 0:
            if args.verbose:
                print "writing 10,000,000 reads"
            for entry in read_dict:
                if args.combine:
                    read_stats[entry] += len(read_dict[entry]) / 4
                    with open(output_dict[entry] + ".fastq", "w") as output:  # file name availability checked before read read-in
                        print>>output, "\n".join(map(str, read_dict[entry]))
                        # print "\t".join(map(str, [output_dict[entry] + ".fastq", len(read_dict[entry]) / 4]))
                    read_dict[entry] = []  # reset read_dict to only accept new reads
                else:
                    for read_dir in read_dict[entry]:
                        read_stats[entry] += len(read_dict[entry][read_dir]) / 4
                        with open(output_dict[entry] + read_dir + ".fastq", "w") as output:  # file name availability checked before read read-in
                            print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
                            # print "\t".join(map(str, [output_dict[entry] + read_dir + ".fastq", len(read_dict[entry]) / 4]))
                        read_dict[entry][read_dir] = []

# need final write for <10million reads at end of file
for entry in read_dict:
    if args.combine:
        read_stats[entry] += len(read_dict[entry]) / 4
        with open(output_dict[entry] + ".fastq", "w") as output:  # file name availability checked before read read-in
            print>>output, "\n".join(map(str, read_dict[entry]))
            # print "\t".join(map(str, [output_dict[entry] + ".fastq", len(read_dict[entry]) / 4]))
        read_dict[entry] = []  # reset read_dict to only accept new reads
    else:
        for read_dir in read_dict[entry]:
            read_stats[entry] += len(read_dict[entry][read_dir]) / 4
            with open(output_dict[entry] + read_dir + ".fastq", "w") as output:  # file name availability checked before read read-in
                print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
                # print "\t".join(map(str, [output_dict[entry] + read_dir + ".fastq", len(read_dict[entry]) / 4]))
            read_dict[entry][read_dir] = []


# Write out new Fastq files
print "Fastq file statistics:"
for entry in read_stats:
    print entry + ".fastq", "\t", read_stats[entry]

# for entry in read_dict:
#     if args.combine:
#         with open(output_dict[entry] + ".fastq", "w") as output:  # file name availability checked before read read-in
#             print>>output, "\n".join(map(str, read_dict[entry]))
#             print "\t".join(map(str, [output_dict[entry] + ".fastq", len(read_dict[entry]) / 4]))
#     else:
#         for read_dir in read_dict[entry]:
#             with open(output_dict[entry] + read_dir + ".fastq", "w") as output:  # file name availability checked before read read-in
#                 print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
#                 print "\t".join(map(str, [output_dict[entry] + read_dir + ".fastq", len(read_dict[entry]) / 4]))
