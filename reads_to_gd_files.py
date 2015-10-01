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


parser = argparse.ArgumentParser(description="read in .fastq.gz file names, generate .gd files. Script allows for single '_' in sample names")
parser.add_argument("-a", "--author", default="Deatherage, Daniel", help="your name")
parser.add_argument("-f", "--fastq", help="direcotry contianing gzipped read files")
parser.add_argument("-r", "--reference", default="genomes/reference/REL606.6.gbk", help="absolute location and name of reference file")
parser.add_argument("-b", "--barrick_backup", action='store_true', help="files are stored in backup.barricklab")
parser.add_argument("-s", "--sra", action='store_true', help="files are on the NCBI SRA archive")
parser.add_argument("-v", "--verbose", action='store_true', help="increase output")
parser.add_argument("-m", "--meta", help="meta data tsv file with headers. 'sample' required, any of the following allowed: 'population', 'time' (meaning generation), 'treatment', 'clone'")
args = parser.parse_args()


allowable_list = ['TIME', 'POPULATION', 'TREATMENT', 'SAMPLE']
meta = {}
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
master_dict = {}
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

sample_names = {}

for files in os.listdir(args.fastq):
    if files.endswith(".fastq.gz"):
        name_parts = files.split("_")
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
        else:
            print "unable to identify sample name for following fastq file: %s" % files
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
    print len(sample_names) + " samples lack meta data. " + sample_names.keys()

#insert somthing here to warn people what gd files you are going to make

loc_of_ref = args.reference
loc_of_reads = args.fastq
if args.barrick_backup:
#    loc_of_ref = loc_of_ref.replace("/home/lab/", "")
#    loc_of_reads = loc_of_reads.replace("/home/lab", "")
    loc_of_ref = loc_of_ref.replace("/corral-repl/utexas/breseq/", "")
    loc_of_reads = loc_of_reads.replace("/corral-repl/utexas/breseq/", "")

for sample in master_dict:
    file_name = "%s.gd" % sample
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
    file_name = "%s.gd" % sample
    with open(file_name, "a") as output:
        print>>output, "#=GENOME_DIFF 1.0"
        print>>output, "#=AUTHOR\t%s" % args.author
        if args.barrick_backup:
            print>>output, "#=REFSEQ\tBarrickLab-Private:%s" % loc_of_ref
            for entry in sample_names[sample]:
                print>>output, "#=READSEQ\tBarrickLab-Private:%s/%s" % (loc_of_reads, entry)
        elif args.sra:
            print>>output, "#=REFSEQ\tBarrickLab-Private:%s" % loc_of_ref  # currently requiring ref to be in private references ... this should be revisited
            for entry in sample_names[sample]:
                print>>output, "#=READSEQ\tSRA:%s" % entry
