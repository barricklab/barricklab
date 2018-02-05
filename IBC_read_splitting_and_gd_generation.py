#!/usr/bin/env python
__author__ = 'ded'
"""
Steps:
    1. Read meta data file
    2. Identify pooled read files from directory
        2a. Unzip files if not already
    3. Split reads based on presence of IBC
    4. Gzip files if they were already
    5. Generate .gd file
        5a. MUST include adaptseq line specific to IBC
        5b. Check if IBC file already exists before generating new
    6. Print recommendations of where to move what files

Meta data file:
Required Columns:
1. Sample
    name of sample, can be descriptive or stain ID (ie DED1), should include at most 1 '_' symbol
2. IBC
    internal barcode sequence used, or IBC1, IBC2 etc for adapter set used
3. Pool
    name of pooled R1 fastq files sample is in, can optionally include .gz, ending 
4. JA_Number
    UTGSAF JA number of pool. Storage of sequencing files uses this information and will be used in the generation of .gd files
5. Ref
    Name of reference file to be used. If stored in barricklab corral, specify as only file name, otherwise website can be specified by beggining with "http" 

Optional Columns:
Meta data to be included in .gd file (ie time, population, treatment, clone)

Expectations:
    1. Reads downloaded from amazon server using _______
    2. Reads moved to Project_###### folder on corral
    3. Script run within corral on TACC
        """

import argparse
import os
import re
import subprocess
from collections import defaultdict
import warnings
import itertools
import datetime
import operator
import sys
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="Read paired end fastq files in, assign reads to different files based on internal barcodes and generate .gd file")
parser.add_argument("-m", "--meta", required=True, help="meta data file containing required columns Sample Name, IBC used, Pooled fastq file, Adapter fasta file")
parser.add_argument("-f", "--fastq", required=True, help="location of pooled fastq files named in meta data file")

# Read splitting options
# parser.add_argument("-r1p", "--read_1_prefix", help="expected sequence in front of IBC on read 1", default="")
# parser.add_argument("-r2p", "--read_2_prefix", help="expected sequence in front of IBC on read 2", default="")
# parser.add_argument("-o", "--offset", help="length of bases to ignore on R1 and R2 after barcode region", default=0, type=int)
parser.add_argument("-c", "--combine", help="combine R1 and R2 into a single fastq output file", default=False, action='store_true')
parser.add_argument("-d", "--distance", help="number of mismatches allowed in R1 or R2 ibc when other ibc is perfect match to expected. default: 0, requires both R1 and R2 ibc match perfectly", default=0)
# parser.add_argument("-p", "--perfect", help="allow assignment of reads to barcodes, only if perfect match", action='store_true')  # depreciated into args.distance

# Genome Diff options
parser.add_argument("-a", "--author", help="Author to be listed in genome diff files", default="IBC_split_and_GD_gen_script")

# Output options
parser.add_argument("-t", "--UnkBCThreshold", help="barcodes exceeding this fraction of the sample with the lowest read count in a given pool will be reported. Set to 1 to report no barcodes; Setting to 0 not advised, will report all barcodes. Lower thresholds will significantly increase output to log file particularly if mismatches are not allowed (distance option set to 0. Default = 0.05", default=0.05, type=float)
parser.add_argument("-v", "--verbose", help="print progress of fastq read in", action='store_true')
parser.add_argument("--log", help="Log file for read statistics, will not catch verbose output", default="barcode_assignment.log.txt")

# other options
parser. add_argument("--testing", help="only work with 50,000 reads for script testing purposes", action='store_true')

# show help if no arguments passed
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# functions
def log(string_to_print, force=False, location=args.log):
    """function to print multiple places at once depending on options"""
    if location:
        with open(location, "a") as print_loc:
            print>>print_loc, string_to_print
    else:
        print string_to_print
    if args.verbose or force:
        print string_to_print


# checks & initiations
assert args.fastq == os.getcwd(), "untested if this works specifying a different location than current."  # TODO test this
assert args.UnkBCThreshold >= 0.0 and args.UnkBCThreshold <= 1, "Specificed unknown barcode threshold: %s; must be between 0 and 1. Lower values will result in more barcodes being reported." % args.UnkBCThreshold
log("\nScript executed at: " + str(datetime.datetime.now()) + "\nWith the following command line options:", force=True)
log("\t" + "\n\t".join(map(str, [str(arg) + "\t" + str(getattr(args, arg)) for arg in vars(args)])), force=True)  # log all input options
if args.testing:
    if not args.verbose:
        args.verbose = True
        log("Verbose feature turned on as required for testing.", force=True)

# Step 1: Read meta data file
with open(args.meta,"r") as meta:
    meta_dict = {}
    pooled_fastqs = []
    first_line = True
    for line in meta:
        line = line.rstrip().split("\t")
        if first_line:
            first_line = False
            header = [y.lower() for y in line]
            assert len(header) >= 5 and header[:5] == ['sample', 'ibc', 'pool', 'ja_number', 'ref'], "Error in meta data header line.\nPlease make sure first 5 columns are ['sample', 'ibc', 'pool', 'ja_number', 'ref']"
            assert len(re.findall("_R1_", line[2])) == 1, "Read 1/2 of pooled reads %s not identifiable as it does not contain'_R1_'. Please verify you have entered the correct pooled fastq file name, and contact Dan if you have" % line[2]
            if len(header) > 5:
                for data in header[5:]:
                    assert data.lower() in ['time', 'population', 'treatment', 'sample', 'clone'], "Error in meta data header line. %s not recognized as meta data.\nPlease make sure optional meta data information only includes:\n'time', 'population', 'treatment', 'sample', or 'clone'.\nContact Dan if you need additional meta data fields" % data
            continue

        assert line[0] not in meta_dict, "%s is given as sample name multiple times in %s.\nWhile there are reasons this may happen, they are very unlikely and this script is not equiped to handle them.\nIf you verify you have identical sample names, contact Dan" %(line[0], args.meta)
        assert len(line) == len(header), "The sample %s has a different number of entries than the header line.\nPlease ensure all samples have all fields filled in.\nIf a sample lacks a value, enter 'None' in that field" % line[0]
        assert re.match("JA", line[3]), "ja_number entry should begin with 'JA' not %s for %s" % (line[3], "\t".join(line))
        meta_dict[line[0]] = dict.fromkeys(header[1:])  # add header keys to new dictionary entry of each sample
        pooled_fastqs = list(set(pooled_fastqs.append(line[2])))  # keep unique list of pooled fastq files

        for number, meta_field in enumerate(header[1:]):
            meta_dict[line[0]][meta_field] = line[number + 1]  # +1 on number so sample isn't included
        assert None not in meta_dict[line[0]].values(), "Header value not found for sample. Coding bug, contact Dan.\n%s" % meta_dict[line[0]]

        #Verify output locations available
        assert not os.path.exists(line[0] + ".gd"), "%s file already exists. Move existing file to new location or delete" % (line[0] + ".gd")
        storage_loc = "/corral-repl/utexas/breseq/genomes/utexas_gsaf/" + "Project_" + line[4]
        if args.combine:
            assert not os.path.exsists(line[0] + ".fastq"), "%s file already exists locally. Move existing file to new location or delete" % (line[0] + ".fastq")
            assert not os.path.exists(storage_loc + line[0] + ".fastq"), "%s file already exists on TACC in the same project folder. This suggests the pool may have already been split. Alternatively, it may mean additional sequencing performed of same samples. Additional sequencing of single samples should not be split in this manner as the generation of .gd files is more complicated. Consult Dan." % (line[0] + ".fastq")
        else:
            for _ in ["_R1", "_R2"]:
                assert not os.path.exists(line[0] + _ + ".fastq"), "%s file already exists locally. Move existing file to new location or delete" % (line[0] + _ + ".fastq")
            assert not os.path.exists(storage_loc + line[0] + _ + ".fastq"), "%s file already exists on TACC in the same project folder. This suggests the pool may have already been split. Alternatively, it may mean additional sequencing performed of same samples. Additional sequencing of single samples should not be split in this manner as the generation of .gd files is more complicated. Consult Dan." % (line[0] + ".fastq")

# Step 2: Verify all listed fastq files in meta data are present in listed directory
fastq_status = {}  #True = gzipped
for fastq_file in pooled_fastqs:
    if fastq_file.endswith('.gz'):
        if fastq_file in os.listdir(args.fastq):
            fastq_status[fastq_file] = True
        if fastq_file[:-3] in os.listdir(args.fastq):
            fastq_status[fastq_file] = False
            log("WARNING: " + fastq_file[:-3] + " already found in " + args.fastq + " and is assumed to be the same file as: " + fastq_file, force=True)
    elif fastq_file + '.gz' in os.listdir(args.fastq):
        fastq_status[fastq_file] = True
        if fastq_file in os.listdir(args.fastq):
            fastq_status[fastq_file] = False
            log("WARNING: " + fastq_file + " already found in " + args.fastq + " and is assumed to be the same file as: " + fastq_file + '.gz', force=True)
    elif fastq_file in os.listdir(args.fastq):
        fastq_status[fastq_file] = False

    assert fastq_file in fastq_status, "Pooled fastq file %s could not be found in %s as a gzipped or non gzipped file" % (fastq_file, args.fastq)

# Step 2a: unzip files if zipped
for fastq_file in fastq_status:
    if fastq_status[fastq_file]:
        if not fastq_file.endswith('.gz'):  # pooled file can be supplied with or without .gz ending, but .gz is always searched for. Only need to unzip if .gz is present and non-.gz isn't
            fastq_file += '.gz'

        #gunzip R1
        gunzip_command = ["gunzip", str(args.fastq) + "/" + str(fastq_file)]
        external_gunzip_call = subprocess.Popen(gunzip_command, stdout=open(args.log,"a"), stderr=open(args.log, "a"))
        external_gunzip_call.wait()

        #gunzip R2
        gunzip_command = ["gunzip", str(args.fastq) + "/" + str(fastq_file).replace("_R1_", "_R2_")]
        external_gunzip_call = subprocess.Popen(gunzip_command, stdout=open(args.log, "a"), stderr=open(args.log, "a"))
        external_gunzip_call.wait()

        for sample in meta_dict:
            if meta_dict[sample]['pool'].endswith(".gz"):
                meta_dict[sample]['pool'] = meta_dict[sample]['pool'][:-3]  # remove .gz ending from file after file expanded
    else:
        pass  # file is already unzipped

# Step 3: Split reads based on IBCs

#get unique list of pooled files to be split (that have all been decompressed appropriately)
all_pools = []
for sample in meta_dict:
    assert not meta_dict[sample]['pool'].endswith(".gz"), "pool appears to still be compressed. Unclear what would cause this, the for loop in 2a block should have changed this"
    assert os.path.exists(args.fastq + "/" + sample), "File: %s\nNot found at: %s" % (meta_dict[sample]['pool'], args.fastq)
    assert os.path.exists(args.fastq + "/" + sample.replace("_R1_", "_R2_")), "Paried end file: %s\nNot found at: %s\nThis script not equipped to handle single end IBC demultiplexing." % (meta_dict[sample]['pool'], args.fastq.replace("_R1_", "_R2_"))

    all_pools.append(meta_dict[sample]['pool'])
all_pools = list(set(all_pools))

# Set up statistic dictionaries
read_counts = {}
mismatched_barcodes = {}
unknown_barcodes = {}

for pool in all_pools:
    unknown = str(pool) + "_unknownBC"
    assert not os.path.exists(unknown + ".fastq"), "unknown barcode assignment file already exists for %s. Move existing file to new location or delete" % pool
    sample_dict = {unknown: unknown}  # reset sample dictionary for each pool. sample_dict[ibc] = samplename for everything besides unknown.

    # Verify pool isn't already in use for statistic dictionaries. This is all very redundant.
    assert pool not in read_counts, "duplicate pool name: %s. This should have been caught earlier in script. Contact Dan" % pool
    assert pool not in mismatched_barcodes, "duplicate pool name: %s. This should have been caught earlier in script. Contact Dan" % pool
    assert pool not in unknown_barcodes, "duplicate pool name: %s. This should have been caught earlier in script. Contact Dan" % pool

    # Set up statistic dictionaries for pool
    read_counts[pool] = defaultdict(int)
    mismatched_barcodes[pool] = 0
    unknown_barcodes[pool] = defaultdict(int)

    # Reset read dictionaries based on desired output. Unknown will be common to all pools. NOTE: unknown is variable based on pool name, NOT STRING "UNKNOWN"
    if args.combine:
        read_dict = {unknown: []}  # generate read dict for all applicable
    else:
        read_dict = {unknown: {"_R1": [], "_R2":[]}}

    # Set up IBCs as keys for all samples expected to be in this pool
    for sample in meta_dict:
        if meta_dict[sample]['pool'] == pool:  #sample is part of pool
            ibc = meta_dict[sample]['ibc'].upper()
            assert ibc not in read_dict, "identical IBC requested for 2 samples in the same pooled sample. If accurate, such samples can not be deconvoluted. Verify accuracy of meta file: %s\n%s has been assigned to multiple samples within pooled fastq file %s" % (args.meta, ibc, pool)
            sample_dict[ibc] = sample
            if args.combine:
                read_dict[ibc] = []
            else:
                read_dict[ibc] = {"_R1": [], "_R2": []}

    # Verify all IBCs are same length. Different length IBCs in single pool very difficult to identify
    assert len(list(set([len(x) for x in read_dict.keys() if x is not unknown]))) == 1, "unequal IBC lengths and can create problems. IBCs are:\n%s\nfor: %s" % (read_dict.keys(), pool)
    bc_length = len(read_dict.keys()[0])

    # Determine distances between barcodes, and if mismatches can be tolerated to level specified at command line
    bc_distances = []
    for pair in itertools.combinations([x for x in read_dict.keys() if x is not unknown], 2):
        bc_distances.append(sum(ch1 != ch2 for ch1, ch2 in zip(pair[0], pair[1])))
        assert bc_distances[-1] > 1, "Distance between IBCs %s and %s is only 1. Single base change can cause misassignment of barcodes. Consult Dan." % (pair[0], pair[1])

    if min(bc_distances) < args.distance:
        new_dist = min(bc_distances)
        message = "2 barcodes are separated by %s, less than user supplied tolerance (%s). Tolerated distance changed to (%s). Other ibc still must match perfectly. WARNING: unsure that this is the best way to do this." % (min(bc_distances), args.distance, new_dist)
        warnings.warn(message)
        args.distance = new_dist

    # Read fastq files
    line_count = 0
    with open(args.fastq + "/" + pool, "r") as fastq1, open(args.fastq + "/" + pool.replace("_R1_", "_R2_"), "r") as fastq2:
        for line1 in fastq1:
            line1 = line1.rstrip()
            line2 = fastq2.next().rstrip()  # advance read 2 file to keep in sync
            line_count += 1
            if line_count % 4 == 1:  # line = header
                header1 = line1
                header2 = line2
                assert (header1.split(" ")[0] == header2.split(" ")[0]) and (len(header1.split(" ")) == len(header1.split(" "))) and len(header1.split(" ")) == 2, "Header lines of R1 and R2 do not agree with expectations: single space in header line, and header before space is identical for both reads.\nHeaders:\n%s\n%s" % (header1, header2)  # this may be inaccurate for different illumina outputs
            elif line_count % 4 == 2:  # line = read
                assert re.match('^[ACTGN]*$', line1) and re.match('^[ACTGN]*$', line2), "Non-sequence characters found on one of th following sequence lines:\n%s\n%s" % (line1, line2)
                bc1 = line1[:bc_length]
                bc2 = line2[:bc_length]

                read1 = line1[bc_length:]
                read2 = line2[bc_length:]

            elif line_count % 4 == 3:  # line = optional sequence descriptor. assume line to use 'standard' "+" notation rather than actual descriptor.
                assert line1 == line2 == "+", "line1 or line2 does not display expected '+' sign on line3 here:\nline1: %s\nline2: %s" % (line1, line2)
            elif line_count % 4 == 0:  # line = quality score. End of single read. Remove quality scores from line, assign read to dictionary based on barcodes.
                quality1 = line1[bc_length:]
                quality2 = line2[bc_length:]

                # Barcode on both reads in agreement, either store in read dictionary, or unknown dictionary
                if bc1 == bc2:
                    if args.combine:
                        try:
                            read_dict[bc1].extend([header1, read1, "+", quality1, header2, read2, "+", quality2])
                        except KeyError:  # barcodes not identified, therefore sore as unknown
                            #  NOTE that the typical "+" symbol on line 3 is changed to what the unidentified barcode was, with first half corresponding to read 1 and second half corresponding to read 2
                            read_dict[unknown].extend([header1, read1, bc1, quality1, header2, read2, bc2, quality2])
                            unknown_barcodes[pool][bc1] += 2  # keep track of how many READS have the barcode, here, bc are identical
                    else:  # R1 and R2 to be kept separate in final output
                        try:
                            read_dict[bc1]["_R1"].extend([header1, read1, "+", quality1])
                            read_dict[bc2]["_R2"].extend([header2, read2, "+", quality2])
                        except KeyError:  # barcodes not specified, therefore sore as unknown
                            #  NOTE that the typical "+" symbol on line 3 is changed to what the unidentified barcode was, with first half corresponding to read 1 and second half corresponding to read 2
                            read_dict[unknown]["_R1"].extend([header1, read1, "+" + bc1, quality1])
                            read_dict[unknown]["_R2"].extend([header2, read2, "+" + bc2, quality2])
                            unknown_barcodes[pool][bc1] += 2  # keep track of how many READS have the barcode. here, bc are identical

                # ibc1 and ibc2 are different
                else:
                    if args.distance > 0:  #mismatches tollerated, identify mismatches and store in read dict if below tolerance, or in unknown if not
                        if bc1 in read_dict or bc2 in read_dict:  # at least one of the barcodes is an expected sequence
                            if bc1 in read_dict and bc2 in read_dict:  # divergent barcodes, store as unknown, and count
                                mismatched_barcodes[pool] += 1
                                read_dict[unknown]["_R1"].extend([header1, read1, "+" + bc1, quality1])
                                read_dict[unknown]["_R2"].extend([header2, read2, "+" + bc2, quality2])
                            else:
                                # Identify which barcode matched to an expected sequence
                                if bc1 in read_dict:
                                    match = bc1
                                else:
                                    match = bc2

                                if sum(ch1 != ch2 for ch1, ch2 in zip(bc1, bc2)) <= args.distance:  # if distance less than or equal to what is tolerated store
                                    if args.combine:
                                        read_dict[match].extend([header1, read1, "+", quality1, header2, read2, "+", quality2])
                                    else:  # R1 and R2 to be kept separate in final output
                                        read_dict[match]["_R1"].extend([header1, read1, "+", quality1])
                                        read_dict[match]["_R2"].extend([header2, read2, "+", quality2])

                                else:  # 2nd barcode more distant than tolerated, both reads stored as unknown
                                    #TODO consider if unknown barcodes should include used IBC that didn't get assigned due to mismatch of 2nd BC. would be a way to tell how many reads might be saved it allowing individual read to assign, or if increased distance tollerance
                                    if bc1 not in read_dict:
                                        unknown_barcodes[pool][bc1] += 1  # keep track of how many READS have the bacrcode
                                    if bc2 not in read_dict:
                                        unknown_barcodes[pool][bc2] += 1  # BC are different here
                                    if args.combine:
                                        read_dict[unknown].extend([header1, read1, "+", quality1, header2, read2, "+", quality2])
                                    else:
                                        read_dict[unknown]["_R1"].extend([header1, read1, "+" + bc1, quality1])
                                        read_dict[unknown]["_R2"].extend([header2, read2, "+" + bc2, quality2])
                        else:  # neither barcode matches expected sequences
                            unknown_barcodes[pool][bc1] += 1  # keep track of how many READS have the barcode
                            unknown_barcodes[pool][bc2] += 1  # BC are different here
                            if args.combine:
                                read_dict[unknown].extend([header1, read1, "+", quality1, header2, read2, "+", quality2])
                            else:
                                read_dict[unknown]["_R1"].extend([header1, read1, "+" + bc1, quality1])
                                read_dict[unknown]["_R2"].extend([header2, read2, "+" + bc2, quality2])
                    else: # distance allowed is 0, therefore can only store as read if BC are equal and expected, which was already dealt with
                        # TODO consider allowing individual half of reads to assign as singleton? ... would also need to go in above. Reason against singletons is they are measure of accuracy of dual end assignment. If you allow singletons, you may as well require only 1 of 2 barcodes req to assign both reads, provided IBCs aren't mismatched (ie IBC1 and IBC2 both represent samples in the pool. Allowing single IBC to assign both reads creates situation where IBC1 may match exactly, IBC2 may have 1 mismatch, but at some frequency, IBC2 is the correct sample with 1 mismatch, and IBC1 is actually 1,2,3,etc mismatches from actual sequence
                        # TODO Again, consider listing used IBC in unknown barcodes to give indication of how many reads would be gained by allowing singletons or increasing distance.
                        if bc1 not in read_dict:
                            unknown_barcodes[pool][bc1] += 1  # keep track of how many READS have the bacrcode
                        if bc2 not in read_dict:
                            unknown_barcodes[pool][bc2] += 1  # BC are different here
                        if args.combine:
                            read_dict[unknown].extend([header1, read1, "+", quality1, header2, read2, "+", quality2])
                        else:
                            read_dict[unknown]["_R1"].extend([header1, read1, "+" + bc1, quality1])
                            read_dict[unknown]["_R2"].extend([header2, read2, "+" + bc2, quality2])
            if args.verbose and line_count % 200000 == 0:
                log(str(line_count / 4) + "read pairs processed")
                if args.testing:
                    break
            if line_count % 40000000 == 0:
                log("Writing 10,000,000 read pairs")
                for entry in read_dict:
                    read_counts[pool][sample_dict[entry]] += len(read_dict[entry]) / 4  # how many reads are being printed for each sample
                    if args.combine:
                        with open(sample_dict[entry] + ".fastq", "a") as output:  # file name availability checked before read read-in
                            print>> output, "\n".join(map(str, read_dict[entry]))
                        read_dict[entry] = []  # reset read_dict to only accept new reads
                    else:
                        for read_dir in read_dict[entry]:
                            with open(sample_dict[entry] + read_dir + ".fastq", "a") as output:  # file name availability checked before read read-in
                                print>> output, "\n".join(map(str, read_dict[entry][read_dir]))
                            read_dict[entry][read_dir] = []  # reset read_dict to only accept new reads


        # Write final/all reads
        log("Writing remaining %i read pairs" % ((line_count % 40000000) / 4))
        for entry in read_dict:
            read_counts[pool][sample_dict[entry]] += len(read_dict[entry]) / 4  # how many reads are being printed for each sample
            if args.combine:
                with open(sample_dict[entry] + ".fastq", "a") as output:  # file name availability checked before read read-in
                    print>> output, "\n".join(map(str, read_dict[entry]))
                read_dict[entry] = []  # reset read_dict to only accept new reads
            else:
                for read_dir in read_dict[entry]:
                    with open(sample_dict[entry] + read_dir + ".fastq", "a") as output:  # file name availability checked before read read-in
                        print>> output, "\n".join(map(str, read_dict[entry][read_dir]))
                    read_dict[entry][read_dir] = []  # reset read_dict to only accept new reads

    log(str(pool) + " read split complete")

# Step 4: gzip split files, and originals
# Better to do within same for loop as the read in, after final fastq write.

    #zip pools
    gzip_command = ["gzip", pool]
    assert not os.path.exists(output_file_name + ".gz"), "re-zipped R1 pool %s already exists" % pool
    external_gzip_call = subprocess.Popen(gzip_command, stdout=open(args.log, "a"), stderr=open(args.log, "a"))
    external_gzip_call.wait()


    gzip_command = ["gzip", pool.replace("_R1_", "_R2_")]
    assert not os.path.exists(output_file_name + ".gz"), "re-zipped R2 pool %s already exists" % pool.replace("_R1_", "_R2_")
    external_gzip_call = subprocess.Popen(gzip_command, stdout=open(args.log, "a"), stderr=open(args.log, "a"))
    external_gzip_call.wait()

    log(str(pool) + " compressed.")

    #zip split files
    for entry in read_dict:
        if args.combine:
            if entry is not "unknown":
                output_file_name = str(sample_dict[entry]) + ".fastq"
            else:
                output_file_name = str(pool) + "_unknown.fastq"

            gzip_command = ["gzip", output_file_name]
            assert not os.path.exists(output_file_name + ".gz"), "re-zipped split_file %s already exists" % output_file_name
            external_gzip_call = subprocess.Popen(gzip_command, stdout=open(args.log, "a"), stderr=open(args.log, "a"))
            external_gzip_call.wait()
            log(output_file_name + " compressed.")

        else:
            for read_dir in read_dict[entry]:
                if entry is not "unknown":
                    output_file_name = str(sample_dict[entry]) + read_dir + ".fastq"
                else:
                    output_file_name = str(pool) + read_dir + "_unknown.fastq"

                gzip_command = ["gzip", output_file_name]
                assert not os.path.exists(output_file_name + ".gz"), "re-zipped split_file %s already exists" % output_file_name
                external_gzip_call = subprocess.Popen(gzip_command, stdout=open(args.log, "a"), stderr=open(args.log, "a"))
                external_gzip_call.wait()
                log(output_file_name + " compressed.")

log("All pooled Fastq files split.\nFastq statistics are as follows:")

# Step 4b: report stats about splitting
Final_Warning = defaultdict(list)
for pool in read_counts:
    log(pool + " had the following samples and read counts:")
    for entry in [z for z in read_counts[pool] if not z.endswith("_unknownBC")] + [k for k in read_counts[pool] if k.endswith("_unknownBC")]:  # make unknown print last for each pool
        log("\t" + entry + "\t" + str(read_counts[pool][entry]))

    # List how many times barcode 1 and 2 suggested assignment to 2 different samples. This is expected to happen at a very low frequency.
    #TODO consider estimating what is a tolerable frequency? (ie mutation rate * hamming distance OR look at several pooled splittings for estimates of actual rates
    log("Number of mismatched read pairs: %s\nThese are read pairs where the internal barcode on read 1 matched a different sample than the internal barcode on read 2.\nIf this number is not very small, contact Dan.\n" % mismatched_barcodes[pool])

    # List barcodes above threshold
    min_read_count = float(min([read_counts[pool][x] for x in read_counts[pool]]))
    if args.UnkBCThreshold < 1:
        log("The following barcodes were identified that do not correspond to an expected sample, and are greater than %s % of the lowest represented sample in this pool." % args.UnkBCThreshold)
        for unkbc in sorted(unknown_barcodes[pool], key=operator.itemgetter(1), reverse=True):
            fraction = unknown_barcodes[pool][unkbc] / min_read_count
            if fraction > args.UnkBCThreshold:
                log("\t" + unkbc + "\t" + unknown_barcodes[pool][unkbc] + "\t" + fraction)
            else:
                # dictionary was sorted in decreasing order, meaning once below threshold, all future will be as well, therefore dont keep iterating
                break
    # Identify barcodes which were more common than sample barcodes
    #TODO consider if this should be done at a given threshold of minimum count (ie count something in the 90th percentile of the minimum read count)
    for unkbc in sorted(unknown_barcodes[pool], key=operator.itemgetter(1), reverse=True):
        if unknown_barcodes[pool][unkbc] > min_read_count:
            log("!!!!! WARNING !!!!!\nThe barcode %s was detected more times than one of the samples from this pool. This suggests major errors such as potentially misentering expected barcode sequence in %s. Verify sequences immediately." % (unkbc, pool), force=True)
            Final_Warning[pool].append(unkbc)  # store unknown IBC in dictionary for reporting at end of script also as this represents a potential big problem with data.
            #TODO consider aborting script here.
        else:
            # dictionary was sorted in decreasing order, meaning one falls below threshold, all future will be as well, therefore dont keep iterating
            break

# Step 5: Generate .gd files & adapter files if new
log("Generating the following genome diff files:")
adapter_files_to_move = []
ref_files_unverified = []
for sample in meta_dict:
    log("\t%s.gd" % sample)
    adapter_file = "/corral-repl/utexas/breseq/genomes/adapters/pooled_sample_IBC_" + meta_dict[sample][ibc].upper() + "_adapter_trim.fa"
    if not os.path.exists(adapter_file) or not os.path.exists(os.getcwd() + "/pooled_sample_IBC_" + meta_dict[sample][ibc].upper() + "_adapter_trim.fa"):  # check if file exists in standard tacc location, or has already been created locally. Second half keeps from making same adapter files for each pool when they are the same (ie same IBC used in multiple pools).
        adapter_files_to_move.append("pooled_sample_IBC_" + meta_dict[sample][ibc].upper() + "_adapter_trim.fa")  # use this to detail where to move what files
        with open("pooled_sample_IBC_" + meta_dict[sample][ibc].upper() + "_adapter_trim.fa", 'W') as adapter_output:
            print>>adapter_output, ">Prefix_Barrick_Internal_Barcode_Pool_" + meta_dict[sample][ibc].upper() + "_Degenerate_External_Barcode/1"
            print>>adapter_output, Seq(meta_dict[sample][ibc].upper() + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG").reverse_complement()
            print>>adapter_output, ">Prefix_Barrick_Internal_Barcode_Pool_" + meta_dict[sample][ibc].upper() + "_Constant_Barcode/2"
            print>>adapter_output, Seq(meta_dict[sample][ibc].upper() + "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT").reverse_complement()
            print>>adapter_output, ">Barrick_Internal_Barcode_Pool_" + meta_dict[sample][ibc].upper() + "_Degenerate_External_Barcode_R1_End/1"
            print>>adapter_output, meta_dict[sample][ibc].upper() + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
            print>>adapter_output, ">Barrick_Internal_Barcode_Pool_" + meta_dict[sample][ibc].upper() + "_Constant_Barcode_R2_End/2"
            print>>adapter_output, meta_dict[sample][ibc].upper() + "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

            # TODO determine if extra A/T needed for library prep, figure out where/when RC of IBC is used and when its used as is

        os.chmod("pooled_sample_IBC_" + meta_dict[sample][ibc].upper() + "_adapter_trim.fa", 0755)  # Modify file permissions to 755, note os.chmod req permissions given in octal hence 0755 to achieve 755

    with open(meta_dict[sample_dict] + ".gd", 'W') as gdoutput:
        print>>gdoutput, "#=GENOME_DIFF 1.0"
        print>>gdoutput, "\t".join(["#=AUTHOR", args.author])

        #Reference
        if re.match("http", meta_dict[sample]['ref']):  # reference is listed on website not corral
            print>>gdoutput, "\t".join(["#REFSEQ", meta_dict[sample]['ref']])
        else:
            print>>gdoutput, "\t".join(["#REFSEQ", "/corral-repl/utexas/breseq/genomes/reference/" + meta_dict[sample]['ref']])
            if not os.path.exists("/corral-repl/utexas/breseq/genomes/reference/" + meta_dict[sample]['ref']):
                ref_files_unverified.append(meta_dict[sample]['ref'])

        #Adapter File
        print>>gdoutput, "\t".join(["#=ADAPTSEQ", adapter_file])

        #Reads
        if args.combine:
            print>>gdoutput, "\t".join(["#=READSEQ", "/corral-repl/utexas/breseq/genomes/Project_" + meta_dict['ja_number'] + "/" + sample + ".fastq.gz"])
        else:
            print>>gdoutput, "\t".join(["#=READSEQ", "/corral-repl/utexas/breseq/genomes/Project_" + meta_dict['ja_number'] + "/" + sample + "_R1.fastq.gz"])
            print>>gdoutput, "\t".join(["#=READSEQ", "/corral-repl/utexas/breseq/genomes/Project_" + meta_dict['ja_number'] + "/" + sample + "_R2.fastq.gz"])

        #Additional Meta Data
        for entry in meta_dict[sample]:
            if entry not in ['sample', 'ibc', 'pool', 'ja_number', 'ref']:
                print>>gdoutput, "\t".join(["#=" + entry, meta_dict[sample][entry]])

    os.chmod(meta_dict[sample_dict] + ".gd", 0755)  #Modify file permissions to 755, note os.chmod req permissions given in octal hence 0755 to achieve 755


# Step 6: Print recomendations of what to do next (where to move files)

if len(adapter_files_to_move) > 0:
    log("The following new adapter file(s) were generated:\n\t" + "\n\t".join(map(str, adapter_files_to_move)))
    log("\n\nThese files should be moved to '/corral-repl/utexas/breseq/genomes/adapters/' using the `mv -i` command. Conflicts were checked for before file generation, and therefore should only exist if problem encountered when trying to access corral when script was run.")

if len(ref_files_unverified) > 0:
    log("The following non-web based reference file(s) was not found at /corral-repl/utexas/breseq/genomes/reference/:\n\t" + "\n\t".join(map(str, ref_files_unverified)))
    log("\n\nCheck if these files are present at '/corral-repl/utexas/breseq/genomes/reference', and add them as needed. Their listing here suggests they are either not present, or you do not have access to the listed directory.")

if Final_Warning:
    log("!!!!! WARNING !!!!!", force=True)
    for pool in Final_Warning:
        for ibc in Final_Warning[pool]:
            log("%s contained more of %s (an unknown barcode) than one of the expected sample barcodes. See above, and verify expected barcodes listed in %s" % (pool, ibc, args.meta), force=True)
    log("!!!!! WARNING !!!!! Very high probability there is a problem with your data", force=True)

