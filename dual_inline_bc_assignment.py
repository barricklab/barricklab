#!/opt/apps/intel15/python/2.7.9/bin/python
__author__ = 'ded'
"""
Script designed to read in dual internal barcode files, and output reads with internal barcodes stripped from reads.
"""

import re
import os
import argparse
from collections import defaultdict
import itertools

# block of improvements to be made:
# TODO: output unknown.fastq changed to reference what the source material was or allow specific name in expected file
# TODO: code combined output fatsqs for mismatch assignment

parser = argparse.ArgumentParser(description="Read paired end fastq files in, assign reads to different files based on barcodes")
parser.add_argument("-l", "--length", default=6, help="length of internal barcodes", type=int)
parser.add_argument("-f1", "--fastq1", help="fastq read1 file to check", required=True)
parser.add_argument("-f2", "--fastq2", help="fastq read2 file to check", required=True)
parser.add_argument("-e", "--expected", help="2 column tsv file, no headers of output file name, and expected barcode sequence. Barcode sequence should be R1R2 no spaces or punctuation.", required=True)
parser.add_argument("-r1p", "--read_1_prefix", help="expected sequence in front of IBC on read 1", default="")
parser.add_argument("-r2p", "--read_2_prefix", help="expected sequence in front of IBC on read 2", default="")
parser.add_argument("-o", "--offset", help="length of bases to ignore on R1 and R2 after barcode region", default=0, type=int)
parser.add_argument("-c", "--combine", help="combine R1 and R2 into a single fastq output file", default=False, action='store_true')
parser.add_argument("-v", "--verbose", help="print progress of fastq read in", action='store_true')
parser.add_argument("-p", "--perfect", help="allow assignment of reads to barcodes, only if perfect match", action='store_true')
args = parser.parse_args()

assert len(args.read_2_prefix) == 0, "sequence provided as read 2 prefix. This is not currently coded as we have no primers which cause this. on expected.tsv file read in, split IBC based on length, insert R2 between 1st and 2nd half. Additional changes to R2 slicing need to be made similar to that of R1 both in dictionary naming and read storing"

# Set up read dictionary
read_dict = {}
output_dict = {"unknown": "unknown_barcodes"}
read_stats = {"unknown": defaultdict(int)}
unknown_barcodes = defaultdict(int)
with open(args.expected, "r") as f:
    for line in f:
        line = line.rstrip().split("\t")
        assert len(line) == 2, "more than 2 columns detected %s" % line
        assert len(line[1]) == 2 * args.length, "barcode given (%s) is different than length (%i) of expected barcode" % (line[1], 2 * args.length)
        if len(args.read_1_prefix) > 0:
            line[1] = args.read_1_prefix + line[1]  # add read_1_prefix in front of listed barcode
        assert line[1] not in read_dict, "identical barcodes given %s" % line[1]
        read_stats[line[1]] = defaultdict(int)
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
            if not args.perfect:
                assert not os.path.exists("junk_sequences" + _ + ".fastq"), "%s file already exists. Move existing file to new location or delete" % (output_dict[entry] + _ + ".fastq")

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
            bc = line1[:args.length + len(args.read_1_prefix)] + line2[:args.length]
            read1 = line1[args.length + len(args.read_1_prefix) + args.offset:]
            read2 = line2[args.length + args.offset:]
        elif line_count % 4 == 3:  # line = optional sequence descriptor. assume line to use 'standard' "+" notation rather than actual descriptor.
            assert line1 == line2 == "+", "line1 or line2 does not display expected '+' sign on line3 here:\nline1: %s\nline2: %s" % (line1, line2)
        elif line_count % 4 == 0:  # line = quality score
            quality1 = line1[args.length + len(args.read_1_prefix) + args.offset:]
            quality2 = line2[args.length + args.offset:]
            if args.combine:
                try:
                    read_dict[bc].extend([header1, read1, "+", quality1, header2, read2, "+", quality2])
                except KeyError:  # barcodes not identified, therefore sore as unknown
                    #  NOTE that the typical "+" symbol on line 3 is changed to what the unidentified barcode was, with first half corresponding to read 1 and second half corresponding to read 2
                    read_dict["unknown"].extend([header1, read1, "+" + bc, quality1, header2, read2, bc, quality2])
                    unknown_barcodes[bc] += 1  # keep track of what barcodes are seen and how many times
            else:  # R1 and R2 to be kept separate in final output
                try:
                    read_dict[bc]["_R1"].extend([header1, read1, "+", quality1])
                    read_dict[bc]["_R2"].extend([header2, read2, "+", quality2])
                except KeyError:  # barcodes not specified, therefore sore as unknown
                    #  NOTE that the typical "+" symbol on line 3 is changed to what the unidentified barcode was, with first half corresponding to read 1 and second half corresponding to read 2
                    read_dict["unknown"]["_R1"].extend([header1, read1, "+" + bc, quality1])
                    read_dict["unknown"]["_R2"].extend([header2, read2, "+" + bc, quality2])
                    unknown_barcodes[bc] += 1  # keep track of what barcodes are seen and how many times
        if args.verbose and line_count % 200000 == 0:
            print line_count / 4, "reads processed"
            # break  # Uncomment for testing subset of reads rather than full read list
        if line_count % 40000000 == 0:
            if args.verbose:
                print "Writing 10,000,000 reads"
            for entry in read_dict:
                if args.combine:
                    read_stats[entry][0] += len(read_dict[entry]) / 4  # 0 because these are prefect matches
                    with open(output_dict[entry] + ".fastq", "a") as output:  # file name availability checked before read read-in
                        print>>output, "\n".join(map(str, read_dict[entry]))
                    read_dict[entry] = []  # reset read_dict to only accept new reads
                else:
                    for read_dir in read_dict[entry]:
                        read_stats[entry][0] += len(read_dict[entry][read_dir]) / 4  # 0 because these are prefect matches
                        with open(output_dict[entry] + read_dir + ".fastq", "a") as output:  # file name availability checked before read read-in
                            print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
                        read_dict[entry][read_dir] = []  # reset read_dict to only accept new reads

# need final write for <10million reads at end of file
if line_count < 40000000:
    print "Writing all %i reads" % (line_count / 4)
else:
    print "Writing remaining %i reads" % (line_count % 40000000) / 4
for entry in read_dict:
    if args.combine:
        read_stats[entry][0] += len(read_dict[entry]) / 4  # 0 because these are prefect matches
        with open(output_dict[entry] + ".fastq", "a") as output:  # file name availability checked before read read-in
            print>>output, "\n".join(map(str, read_dict[entry]))
        read_dict[entry] = []  # reset read_dict to only accept new reads (important if mismatches to be allowed)
    else:
        for read_dir in read_dict[entry]:
            read_stats[entry][0] += len(read_dict[entry][read_dir]) / 4  # 0 because these are prefect matches
            with open(output_dict[entry] + read_dir + ".fastq", "a") as output:  # file name availability checked before read read-in
                print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
            read_dict[entry][read_dir] = []  # reset read_dict to only accept new reads (important if mismatches to be allowed)

# Write out new Fastq file statistics
print "Fastq file statistics for reads with perfect matches:"
for entry in read_stats:
    if entry == "unknown":
        continue
    print output_dict[entry] + ".fastq", "\t", read_stats[entry][0]  # 0 because only perfect matches have been and will be stored
print output_dict["unknown"] + ".fastq", "\t", read_stats["unknown"][0]  # 0 because only perfect matches have been and will be stored

if not args.perfect:  # mismatches allowed, therefore try to assign reads in the "unknown" group to barcode

    print "\nAttempting assignment of %i unknown barcodes with mistmatches." % read_stats["unknown"][0]

    # Distances to be calculated using Levenshtein Distance rather than hamming distance to allow for indels.
    def levenshtein(s1, s2):
        """
        Calculates levenshtein distance between s1 and s2. see https://en.wikipedia.org/wiki/Levenshtein_distance for more info.
        :param s1: string 1
        :param s2: string 2
        :return: levenshtein distance between two sequences
        variation of code found at https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
        """
        # Coded under assumption that as barcode regions are trimmed from raw reads and are appropriate length.
        # see above links for use with unequal lengths

        previous_row = range(len(s2) + 1)
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                #  TODO evaluate concerns that if BC contains insertions, the remaining read contains adapter sequence to be used for mapping, and figure out how to modify this
                insertions = previous_row[j + 1] + 1  # j+1 instead of j since previous_row and current_row are one character longer than s2
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)  # True has value of 1 False has value of 0
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row

        return previous_row[-1]

    # Determine distance barcodes are from each other
    barcode_distances = []
    for pair in itertools.combinations([x for x in read_dict.keys() if x is not "unknown"], 2):
        barcode_distances.append(levenshtein(pair[0], pair[1]))
    distance_tolerated = (min(barcode_distances) / 2.0) - 1  # observed BC must be more than twice as close to 1 barcode as any other
    assert distance_tolerated >= 1, "distance tollerated <1, therefore, no mutations are allowed ... reads "
    print "All barcodes are at least %i mutations away from one another. Therefore an unknown barcode must be within %i mutations of a known barcode to be assigned. Ties are not possible" % (min(barcode_distances), distance_tolerated)

    # group unknown barcodes with known barcodes if they within given distance
    mismatch_dict = {}
    for unknown_barcode in unknown_barcodes:
        best = float("inf")
        for bc in [x for x in read_dict.keys() if x is not "unknown"]:
            distance = levenshtein(unknown_barcode, bc)
            if distance < best:
                best = distance
                mismatch_dict[unknown_barcode] = bc
            elif distance == best:
                mismatch_dict[unknown_barcode] = [mismatch_dict[unknown_barcode], bc]
        if best > distance_tolerated:
            mismatch_dict[unknown_barcode] = "unknown"
            read_stats["unknown"][distance] += unknown_barcodes[unknown_barcode]  # add number of reads observed with this barcode
            continue
        assert isinstance(mismatch_dict[unknown_barcode], str), "Barcode tied with a distance tolerated for assignment!\nunknown barcode:\t%s\nmatches:\t%s\ndistance:\t%i" % (unknown_barcode, mismatch_dict[unknown_barcode], best)
        read_stats[mismatch_dict[unknown_barcodes[unknown_barcode]]][best] += unknown_barcodes[unknown_barcode]  # add number of reads observed with this barcode based on distance

    # Read unknown barcode.fastq files in.
    if args.combine:
        assert False, "this still needs to be coded. Let Dan know, and in the mean time, repeat analysis without -c option and cat files together after run finishes"
        # this may need to be reading in in sets of 8, or reading in single file and just outputing as it goes or at intervals
        # if combine before, combine now, so sets of 8 not likely needed or beneficial
    else:
        # Fastq files passed assertions prior to writing, therefore limited use here for speed
        line_count = 0
        with open("unknown_barcodes_R1.fastq", "r") as fastq1, open("unknown_barcodes_R2.fastq", "r") as fastq2:  # TODO revisit names if unknown to represent raw input names
            for line1 in fastq1:
                line1 = line1.rstrip()
                line2 = fastq2.next().rstrip()  # advance read 2 file to keep in sync
                line_count += 1
                if line_count % 4 == 1:  # line = header
                    header1 = line1
                    header2 = line2
                elif line_count % 4 == 2:  # line = read
                    read1 = line1
                    read2 = line2
                elif line_count % 4 == 3:  # line = bc
                    assert line1 == line2 and line1.startswith("+") and len(line1) == (args.length * 2 + len(args.read_1_prefix)), "line 3 expected to be same for both reads, expected to start with + and barcode region/length.\nRead1:\t%s\nRead2:\t%s\nLine:\t" % (line1, line2, line_count)
                    bc = line1.lstrip("+")
                elif line_count % 4 == 0:  # line = quality score
                    quality1 = line1
                    quality2 = line2
                    new_bc = mismatch_dict[bc]
                    if new_bc is not "unknown":
                        read_dict[mismatch_dict[bc]]["_R1"].extend([header1, read1, "+", quality1])
                        read_dict[mismatch_dict[bc]]["_R2"].extend([header2, read2, "+", quality2])
                    else:
                        read_dict["unknown"]["_R1"].extend([header1, read1, "+" + bc, quality1])
                        read_dict["unknown"]["_R2"].extend([header2, read2, "+" + bc, quality2])
                if args.verbose and line_count % 200000 == 0:
                    print line_count / 4, "unknown reads processed"
                    # break  # Uncomment for testing subset of reads rather than full read list

                if line_count % 40000000 == 0:
                    if args.verbose:
                        print "Writing 10,000,000 unknown reads"
                    for entry in read_dict:
                        if entry is not "unknown":
                            for read_dir in read_dict[entry]:
                                # TODO stats on new assignments #read_stats[entry] += len(read_dict[entry][read_dir]) / 4
                                with open(output_dict[entry] + read_dir + ".fastq", "a") as output:  # file name availability checked before read read-in
                                    print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
                                read_dict[entry][read_dir] = []  # reset read_dict to only accept new reads
                        else:
                            for read_dir in read_dict["unknown"]:
                                with open("junk_sequences" + read_dir + ".fastq", "a") as output:
                                    print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
                                read_dict["unknown"][read_dir] = []  # reset read_dict to only accept new reads

    # need final write for <10million unknown reads at end of file
    if line_count < 40000000:
        print "Writing all %i unknown reads" % line_count / 4
    else:
        print "Writing remaining %i unknown reads" % (line_count % 40000000) / 4
    for entry in read_dict:
        if entry is not "unknown":
            for read_dir in read_dict[entry]:
                # TODO stats on new assignments #read_stats[entry] += len(read_dict[entry][read_dir]) / 4
                with open(output_dict[entry] + read_dir + ".fastq", "a") as output:  # file name availability checked before read read-in
                    print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
                read_dict[entry][read_dir] = []  # reset read_dict to only accept new reads should be irrelevant
        else:
            for read_dir in read_dict["unknown"]:
                with open("junk_sequences" + read_dir + ".fastq", "a") as output:
                    print>>output, "\n".join(map(str, read_dict[entry][read_dir]))
                read_dict["unknown"][read_dir] = []  # reset read_dict to only accept new reads should be irrelevant

    all_distances = list(set((itertools.chain.from_iterable([read_stats[x].keys() for x in read_stats]))))  # generate list of all "best" distances detected
    all_distances = all_distances.sort()  # if sort included with above, returns no values
    print "\t".join(map(str, ['Sample/Distance'] + all_distances))
    for bc in read_stats:
        to_print = [bc]
        for dist in all_distances:
            try:
                to_print.append(read_stats[bc][dist])
            except KeyError:
                to_print.append(0)
        print "\t".join(map(str, to_print))


print "Script Completed Sucessfully"

# print "5-2-2016"
# print "Based on Brian's tnSeq data, R1 barcode has an additional A as the first base off. Script currently treats this as generalized fact."
# print "If output looks odd (ie 100% of reads going to unk), check this first. Dan and Sean"
# above removed 5-31-2016 based on Dacia's data. For Brian's Data, script now requires -read_1_prefix "A" for brian's transposon library


