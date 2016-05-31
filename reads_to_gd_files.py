#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:52:12 2013
Find names of reads in directory, create gd files based on those names

improvements needed:
    print .gd files in a different location than the read files
    print what names of reads will be, possibly with confirmation for if this is correct
    need ncbi and other support

DO SOMETHING WITH TRAILING / marks on autocomplete with tab for fastq directory
@author: ded
"""

import os
import re
import argparse
import sys
import re

parser = argparse.ArgumentParser(description="read in .fastq.gz file names, generate .gd files. Script allows for single '_' in sample names")
parser.add_argument("-a", "--author", help="Your name.")
parser.add_argument("-f", "--fastq", help="Absolute path to directory containing gzipped read files. Example: /corral-repl/utexas/breseq/genomes/utexas_gsaf/Project_###", required=True)
parser.add_argument("-i", "--index", action='store_true', help="Index files exist in fastq directory, but should be ignored.")
parser.add_argument("-r", "--reference", nargs='*', help="Absolute location and name of reference file. Example: /corral-repl/utexas/breseq/genomes/reference/REL606.6.gbk", default="None")
parser.add_argument("-b", "--barricklab", action='store_true', help="Files are stored in barricklab location on corral.")
parser.add_argument("-s", "--sra", action='store_true', help="Files are on the NCBI SRA archive.")
parser.add_argument("-o", "--output", default="new_gd_files", help="Output directory for .gd files (will be created if it doesn't exist).")
parser.add_argument("-v", "--verbose", action='store_true', help="Increase output")
parser.add_argument("-m", "--meta", help="Optional meta data tsv file with headers. 'sample' required, any of the following allowed: 'population', 'time' (meaning generation), 'treatment', 'clone'")


# show help if no arguments passed
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Check that barricklab or sra is location of read files
if args.barricklab == args.sra == False:
    parser.print_help()
    print "\n\n!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!\nNeither -barricklab nor -sra commands given. 1 must be supplied. Script aborted.\n\n"
    sys.exit(1)
if args.barricklab is not False and args.sra is not False:
    parser.print_help()
    print "\n\n!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!\nBoth -barricklab and -sra commands given. Only 1 may be supplied. Script aborted.\n\n"
    sys.exit(1)

# Check that absolute paths for fastq and reference files are used.
root_warn = [re.search("^/", args.fastq)] + [re.search("None|^/", ref_warn_check) for ref_warn_check in args.reference]
if root_warn.count(None):
    parser.print_help()
    print "\n\n!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!\nPlease specify fastq directory and/or reference file from root.\n"
    print "\n".join(map(str, [args.fastq] + args.reference))
    sys.exit(1)

allowable_list = ['TIME', 'POPULATION', 'TREATMENT', 'SAMPLE', 'CLONE']
meta = {}
sample_names = {}
master_dict = {}

if args.meta is not None:
    
    with open(args.meta, "r") as F:
        for line in F:
            line = line.rstrip().split("\t")
            if len(meta) == 0:
                meta['id'] = [x.upper() for x in line]
                assert 'SAMPLE' in meta['id'], '"sample" name missing from meta file read in'
                for thing in meta['id']:
                    assert thing in allowable_list, '%s is unknown meta data, consider changing name to one of %s if alternative name, or check spelling' % (thing, allowable_list)
                meta['samples'] = []
                continue
            meta['samples'].append(line)

    sample_list = []
    for sample in meta['samples']:
        sample_id = ''
        other = []
        for index, name in enumerate(meta['id']):
            if meta['id'][index] == "SAMPLE":
                sample_id = sample[index]
            else:
                other.append([meta['id'][index], sample[index]])
        master_dict[sample_id] = other

#for sample in master_dict:
#    for thing in master_dict[sample]:
#        print "#=%s\t%s" % (thing[0], thing[1])

for files in os.listdir(args.fastq):
    if files.endswith(".fastq.gz"):
        if args.index:
            if re.search("L00\d_I\d_00\d.fastq.gz", files):
                continue  # ignore fastq files that are generated from index reads.

        sample_name = files
        if re.search('_S\d+_L\d\d\d_R\d_\d\d\d.fastq.gz$', sample_name):  # @DED 3-30-16 new fastq naming system being used by gsaf
            sample_name = re.sub(r'_S\d+_L\d\d\d_R\d_\d\d\d.fastq.gz$', r'', sample_name)
            sample_name = re.sub(r'^JA\d+_', r'', sample_name)  # remove common unnecessary JA IDs
            sample_name = re.sub(r'^Sample_', r'', sample_name)  # remove common Sample_ from file names
            print files, "becomes", sample_name

        else:
            sample_name = sample_name.replace(".fastq.gz", "")
            # ditch things like this _L002_R2_001 ... @DED 3-30-16 this assumes old gsaf file naming system which had 6 letter bar code prior to lane
            sample_name = re.sub(r'_L\d\d\d_R\d_\d\d\d$', r'', sample_name)
            sample_name = re.sub(r'^JA\d+_', r'', sample_name)  # remove common unnecessary JA IDs
            sample_name = re.sub(r'^Sample_', r'', sample_name)  # remove common Sample_ from file names
            print sample_name
        
        name_parts = sample_name.split("_")
            
        if re.match("[ACTG]{6}", name_parts[1]):
            if name_parts[0] in master_dict:
                master_dict[name_parts[0]].append(['READSEQ', files])
            elif name_parts[0] in sample_names:
                sample_names[name_parts[0]].append(files)
            elif name_parts[0] not in sample_names:
                sample_names[name_parts[0]] = [files]
        elif re.match("[ACTG]{6}", name_parts[2]):  # allows for 1 underscore character in name system
            if name_parts[0] + "_" + name_parts[1] in master_dict:
                master_dict[name_parts[0] + "_" + name_parts[1]].append(['READSEQ', files])
            elif name_parts[0] + "_" + name_parts[1] in sample_names:
                sample_names[name_parts[0] + "_" + name_parts[1]].append(files)
            elif name_parts[0] + "_" + name_parts[1] not in sample_names:
                sample_names[name_parts[0] + "_" + name_parts[1]] = [files]
        else:  # @DED 3-30-16 this is what is catching the new gsaf fastq output names
            if sample_name in master_dict:
                master_dict[sample_name].append(['READSEQ', files])
            elif sample_name in sample_names:
                sample_names[sample_name].append(files)
            elif sample_name not in sample_names:
                sample_names[sample_name] = [files]
            
    if files.endswith(".sra") and args.sra:
        with open(files, "r") as F:
            for line in F:
                line = line.rstrip().split(",")
                if line[0] == "Run":
                    continue
                else:
                    x = "%s_1" % line[0]
                    y = "%s_2" % line[0]
                    sample_names[line[27]] = [x, y]
                    #try:
                    #    sample_names[line[27]].append(line[0])
                    #except(KeyError):



if len(sample_names) > 0:
    print str(len(sample_names)) + " samples lack meta data. " + str(sample_names.keys())



if not os.path.exists(args.output):
    os.makedirs(args.output)

loc_of_ref = args.reference
loc_of_reads = args.fastq
if args.barricklab:
#    loc_of_ref = loc_of_ref.replace("/home/lab/", "")
#    loc_of_reads = loc_of_reads.replace("/home/lab", "")
    loc_of_ref = loc_of_ref.replace("/corral-repl/utexas/breseq/", "")
    loc_of_reads = loc_of_reads.replace("/corral-repl/utexas/breseq/", "")


print "The following files are being written:"  # insert somthing here to warn people what gd files you are going to make

for sample in master_dict:
    file_name = "%s/%s.gd" % (args.output.rstrip("/"), sample)
    print "\t", file_name
    with open(file_name, "a") as output:
        print>>output, "#=GENOME_DIFF 1.0"
        print>>output, "#=AUTHOR\t%s" % args.author
        for thing in master_dict[sample]:
            if thing[0] == 'READSEQ':
                print>>output, "#=READSEQ\tBarrickLab-Private:%s/%s" % (loc_of_reads, thing[1])
            else:
                print>>output, "#=%s\t%s" % (thing[0], thing[1])
        print>>output, "#=REFSEQ\tBarrickLab-Private:%s" % loc_of_ref

for sample in sample_names:
    file_name = "%s/%s.gd" % (args.output.rstrip("/"), sample)
    print "\t", file_name
    with open(file_name, "a") as output:
        print>>output, "#=GENOME_DIFF 1.0"
        print>>output, "#=AUTHOR\t%s" % args.author
        if args.barricklab:
            print>>output, "#=REFSEQ\tBarrickLab-Private:%s" % loc_of_ref
            for entry in sample_names[sample]:
                print>>output, "#=READSEQ\tBarrickLab-Private:%s/%s" % (loc_of_reads, entry)
        elif args.sra:
            print>>output, "#=REFSEQ\tBarrickLab-Private:%s" % loc_of_ref  # currently requiring ref to be in private references ... this should be revisited
            for entry in sample_names[sample]:
                print>>output, "#=READSEQ\tSRA:%s" % entry
