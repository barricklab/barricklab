#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:52:06 2013
Duplex Seq Input To Breseq
Step 1
    Identify Molecular Index (MI) on each read pair
Step 2
    Output reads to 1 of 3 intermediate files: Double_MI.fastq, Single_MI.fastq, Zero_MI_pair.fastq. To be based on how many correct MI are identified in the pair.
    MI should be removed from read sequence string, and appended to header
    Quality score should be removed from MI regions
****Future Improvment:
        be more liberal with MI and constant region (ie 11bp then 4bp constant for synthesis error). Need to let non-appropriate reads drive this.
        Quality score of MI region should be used to calculate MI
Step 3
    Focus only on Double_MI.fastq first
    Store tuples in dictionary as


#Breseq does not make mutation calls based on "N" sequences. Therefore it is not necessary to identify 1st round PCR errors. #### NOT TRUE when coverag uneven


would be nice to add something for reading in gd files ot generate command file for better merge with breseq

current stuff printing to screen should default to printing to log file (specifie from argrse)
function use of main()
deletion of intermediate file dictioaries afterwrite
append stats of run (ie histogram of reads per MI, total reduction, number of Ns, etc) possibly multiple files,


if multiple reads from a single sample:
    for to_combine in *L007*.fastq;do name=$(echo $to_combine|sed 's/_.*//');files=$(ls $name*);files=$(echo $files|sed 's/\n/ /'); cat $files > $name.combined.fastq;done
    #above is not tested very well, but header assertions exist to catch any cat errors


for R1_name in Raw_R1/*.fastq;do name=$(echo $R1_name|sed 's/Raw_R1\///'|sed 's/_.*//');R2_name=$(echo Raw_R2/$name*);echo "superread.py -f1 $R1_name -f2 $R2_name -p $name --log $name.log.txt";done >commands

@author: ded
"""

import re
import argparse
import sys
import datetime
import numpy
from collections import defaultdict

parser = argparse.ArgumentParser(description="calculate how many reads have expected duplex seq tag")
parser.add_argument("-c", "--constant_region", default= "CAGT", help="expected constant region sequece")
parser.add_argument("-l", "--length", default = int(12), help="length of random MI region")
parser.add_argument("-f1", "--fastq1", help="fastq read1 file to check")
parser.add_argument("-f2", "--fastq2", help="fastq read2 file to check")
parser.add_argument("-cf", "--cutoff_frequency", default = float(0.5), type=float, help="minimum cutoff frequency required to call concensus base")
parser.add_argument("-p", "--prefix", help = "prefix for output files")
parser.add_argument("-m", "--minimum_reads", default = 2, type=int, help = "minimum number of reads needed to support SSCS reads")
parser.add_argument("-v", "--verbose", action='store_true', help="increase output (not set to do anything)")
parser.add_argument("-i", "--intermediate_fastqs", action='store_true', help="remove MI and output 'double, single, zero'.fastq files based on presence of 2,1,0 MI. DOES NOT EFFECT DCS, SSCS, or Singleton output files.")
parser.add_argument("-d", "--DCS", action='store_true', default = False, help="calculate DCS sequence, off by default. Turns SSCS on as well.")
parser.add_argument("-s", "--SSCS", action='store_true', default = False, help="calculate SSCS sequence, off by default. IF DCS specificed, automatically on")
#parser.add_argument("-t", "--trimmed", help="enter 2 for requiring correctly identified MI on both reads, enter 1 for requiring correct MI on single read. non")#hasto be on for anything to happen?
parser.add_argument("--log", help = "name of output log file")
args = parser.parse_args()

start = datetime.datetime.now()

if args.log:
    log = open(args.log, "w")
else:
    log = sys.stdout

print>>log, "processing sample %s with read files:\n%s\n%s" %(args.prefix, args.fastq1, args.fastq2)
print>>log, "start time: %s" %datetime.datetime.now()
molecular_index_search_string = "^[ACTGN]{%d}%s" %(args.length, args.constant_region)
fastq_dict = {"Double_MI.fastq": [], "Single_MI.fastq": [], "Zero_MI.fastq":[]}
#groups_dict = {"DCS" : {}, "SSCS-R1" : {}, "SSCS-R2" : {}} #depreciated?
#final_output_dict = {"DCS" : {}, "SSCS-1" : {}, "SSCS-2": {}} #depreciated?
master_dict = {}
good_reads = 0
Base_discrepancies = {}
Unresolved_consensus = defaultdict(int)
total_concensus_attempts = 0
base_positions = []



def de_duplex_read(header1, read1, quality1, header2, read2, quality2):
    """Function designed to take in fastq information and output reads lacking molecular index"""
    global good_reads
    if args.intermediate_fastqs:
        if re.match(molecular_index_search_string, read1) and re.match(molecular_index_search_string, read2):
            header1 = header1+"#MI:"+read1[:12]
            read1 = read1[17:]
            quality1 = quality1[17:]

            header2 = header2+"#MI:"+read2[:12]
            read2 = read2[17:]
            quality2 = quality2[17:]

            good_reads+=2
            group_assignment(header1, read1, header2, read2) #quality1, quality2 #moved outside of function as quality ends up too high to use

            fastq_dict["Double_MI.fastq"].extend([[header1, read1, "+", quality1],[header2, read2, "+", quality2]])

        elif re.match(molecular_index_search_string, read1):
            header1 = header1+"#MI:"+read1[:12]
            read1 = read1[17:]
            quality1 = quality1[17:]
            fastq_dict["Single_MI.fastq"].extend([[header1, read1, "+", quality1],[header2, read2, "+", quality2]])

        elif re.match(molecular_index_search_string, read2):
            header2 = header2+"#MI:"+read2[:12]
            read2 = read2[17:]
            quality2 = quality2[17:]
            fastq_dict["Single_MI.fastq"].extend([[header1, read1, "+", quality1],[header2, read2, "+", quality2]])

        else:
            fastq_dict["Zero_MI.fastq"].extend([[header1, read1, "+", quality1],[header2, read2, "+", quality2]])

    else:
        if re.match(molecular_index_search_string, read1) and re.match(molecular_index_search_string, read2):
            header1 = header1+"#MI:"+read1[:12]
            read1 = read1[17:]
            quality1 = quality1[17:]

            header2 = header2+"#MI:"+read2[:12]
            read2 = read2[17:]
            quality2 = quality2[17:]

            good_reads+=2
            group_assignment(header1, read1, header2, read2) #quality1, quality2 #moved outside of function as quality ends up too high to use

def group_assignment(header1, read1, header2, read2): #quality1, quality2 #moved outside of function as quality ends up too high to use
    """function designed to take in fastq information where MI has been added to header, and group reads with identical MI for SingleStrandConensusGroups
    Currently doing nothing with full header or quality score, thisshouldprobably be upgraded to do both things"""

    R1_tag = header1.split('#MI:')[-1]+header2.split('#MI:')[-1]
    R2_tag = header2.split('#MI:')[-1]+header1.split('#MI:')[-1]

    for tag in [R1_tag, R2_tag]:
        if tag not in master_dict:
            master_dict[tag] = {}
            master_dict[tag]["Read1"] = []
            master_dict[tag]["Read2"] = []
            master_dict[tag]["SSCS1"] = []
            master_dict[tag]["SSCS2"] = []
            master_dict[tag]["DCS"] = []
    master_dict[R1_tag]["Read1"].append(read1)
    master_dict[R2_tag]["Read2"].append(read2)


def consensus_generation(list_of_sequences, cutoff, read_length):
    """take in list of sequences, and consensus sequence abve threshold.
    Largely copied from Scott Kennedy's consensusMaker.py script used by Schmitt et all
    partially copied from Scott Kennedy's script
    ok, not really copied from tht script anymore"""
    global total_concensus_attempts
    total_concensus_attempts +=1
    global base_positions
    consensus_read = []
    if len(list_of_sequences) not in Base_discrepancies:
        Base_discrepancies[len(list_of_sequences)] = {}

    for i in xrange(read_length) : #Count the types of nucleotides at a position in a read. i is the nucleotide index within a read in list_of_sequences
        consensus_dict = {'A':0, 'C':0, 'T':0, 'G':0, 'N':0} #set dictionary  to empty for each base within the readgroup
        if i not in Base_discrepancies[len(list_of_sequences)]:
            Base_discrepancies[len(list_of_sequences)][i] = []
        if i not in base_positions:
            base_positions.append(i)
        for j in xrange(len(list_of_sequences)) : #Do this for every read that comprises a SMI group. j is the read index within list_of_sequences
            consensus_dict[list_of_sequences[j][i]] +=1.0

        #add an assertion for the sum total of base counts equal to the number of sequences?
        if consensus_dict[max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key]))]/float(len(list_of_sequences)) >cutoff: #this would fall apart if a cutoff of less than 0.5 is used
            consensus_read.append(max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key])))
            Base_discrepancies[len(list_of_sequences)][i].append(consensus_dict[max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key]))]/float(len(list_of_sequences)))
            if Base_discrepancies[len(list_of_sequences)][i][-1] != 1.0 and len(list_of_sequences) ==2:
                print len(list_of_sequences), Base_discrepancies[len(list_of_sequences)][i][-1]
#            try:
#                Base_discrepancies[i].append(consensus_dict[max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key]))]/float(len(list_of_sequences)))
#            except (KeyError):
#                Base_discrepancies[i] = [consensus_dict[max(consensus_dict.iterkeys(), key=(lambda key: consensus_dict[key]))]/float(len(list_of_sequences))]
        else:
            consensus_read.append("N")
            Unresolved_consensus[i] += 1
    consensus_read = "".join(consensus_read)
    return consensus_read

#read fastq1 an fastq2 in together
line_count = 0
with open(args.fastq1, "r") as Fastq1, open(args.fastq2, "r") as Fastq2:
    for line1 in Fastq1:
        line1 = line1.rstrip()
        line2 = Fastq2.next().rstrip()
        line_count +=1
        if line_count%4 ==1: #line = header
            header1 = line1
            header2 = line2
            assert (header1.split(" ")[0] == header2.split(" ")[0]) and (len(header1.split(" ")) == len(header1.split(" "))) and len(header1.split(" ")) == 2, "Header lines of R1 and R2 do not agree with expectations: single space in header line, and header before space is identical for both reads.\nHeaders:\n%s\n%s" %(header1, header2) #this may be inaccurate for different illumina outputs
        elif line_count%4 ==2:#line = read
            read1 = line1
            read2 = line2
            assert re.match('^[ACTGN]*$', read1) and re.match('^[ACTGN]*$', read2), "Non-sequence characters found on one of th following sequence lines:\n%s\n%s" %(line1,line2)
        elif line_count%4 ==3: #line = optional sequence descriptor. assume line to use 'standard' "+" notation rather than actual descriptor.
            assert line1 == line2 == "+", "line1 or line2 does not display expected '+' sign on line3 here:\nline1: %s\nline2: %s" %(line1, line2)
        elif line_count%4 ==0:
            quality1 = line1
            quality2 = line2
            de_duplex_read(header1, read1, quality1, header2, read2, quality2)
        if line_count%200000 == 0:
            print>>log, line_count/4, "reads processed"
           # break #for testing subsetof reads rather than full read list
sys.stdout.flush() #dump processed reads out before continuing

##### write 3 fastq files from fastq_dict serves as early intermediate of where to focus on future improvements
if args.intermediate_fastqs:
    for entry in fastq_dict:
        output_name = args.prefix+entry
        with open(output_name,'w') as output:
            for read in fastq_dict[entry]:
                print>>output, "\n".join(read)

#generte SSCS and DCS consensus sequences
Reads_per_MI_Dict = defaultdict(int) #usedefault dict to simplif next 6 lines
for MI in master_dict:
#    if len(master_dict[MI]["Read1"]) not in Reads_per_MI_Dict:
#        Reads_per_MI_Dict[len(master_dict[MI]["Read1"])] = 0
#    if len(master_dict[MI]["Read2"]) not in Reads_per_MI_Dict:
#        Reads_per_MI_Dict[len(master_dict[MI]["Read2"])] = 0
    Reads_per_MI_Dict[len(master_dict[MI]["Read1"])] +=1
    Reads_per_MI_Dict[len(master_dict[MI]["Read2"])] +=1

    if args.SSCS or args.DCS:
        if len(master_dict[MI]["Read1"]) > (args.minimum_reads-1):
            assert len(master_dict[MI]["SSCS1"]) == 0, "duplicated MI: %s" %MI
            master_dict[MI]["SSCS1"] = [consensus_generation(master_dict[MI]["Read1"], args.cutoff_frequency, len(master_dict[MI]["Read1"][0]))]
        if len(master_dict[MI]["Read2"]) > (args.minimum_reads-1):
            assert len(master_dict[MI]["SSCS2"]) == 0, "duplicated MI: %s" %MI
            master_dict[MI]["SSCS2"] = [consensus_generation(master_dict[MI]["Read2"], args.cutoff_frequency, len(master_dict[MI]["Read2"][0]))]
        if args.DCS:
            if len(master_dict[MI]["SSCS1"]) ==1 and len(master_dict[MI]["SSCS2"]) == 1:
                assert len(master_dict[MI]["DCS"]) == 0, "duplicated MI: %s" %MI
                master_dict[MI]["DCS"] = [consensus_generation([master_dict[MI]["SSCS1"][0] , master_dict[MI]["SSCS2"][0] ], 0.5, len(master_dict[MI]["Read2"][0]))]

#write DCS, SSCS, singleton sequences
double_singletons = 0
with open(args.prefix+"_DCS.fastq", "w") as DCS, open(args.prefix+"_SSCS.fastq", "w") as SSCS, open(args.prefix+"_singletons.fastq", "w") as singles:#, open(args.prefix+"_Double_Singletons.fastq", "w") as doubles:
    dcs_count = sscs_count = single_count = 0
    for MI in master_dict:
        if len(master_dict[MI]["DCS"]) ==1:
            print>>DCS, "\n".join(map(str, ["@"+str(MI)+"-"+str(dcs_count), master_dict[MI]["DCS"][0], "+", 'I'*len(master_dict[MI]["DCS"][0])]))
            dcs_count +=1
        elif len(master_dict[MI]["SSCS1"]) == 1:#if sscs1 and 2 were both 1 then there ould be a DCS
            print>>SSCS, "\n".join(map(str, ["@"+str(MI)+"-"+str(sscs_count), master_dict[MI]["SSCS1"][0], "+", '5'*len(master_dict[MI]["SSCS1"][0])]))
            sscs_count +=1
        elif len(master_dict[MI]["SSCS2"]) == 1: #if sscs1 and 2 were both 1 then there ould be a DCS
            print>>SSCS, "\n".join(map(str, ["@"+str(MI)+"-"+str(sscs_count), master_dict[MI]["SSCS2"][0], "+", '5'*len(master_dict[MI]["SSCS2"][0])])) #5 = qualit score of 20
            sscs_count +=1
        if len(master_dict[MI]["Read1"]) == 1: # if read1 only has1 read then no sscs or dcs read can exist
            print>>singles, "\n".join(map(str, ["@"+str(MI)+"-"+str(single_count), master_dict[MI]["Read1"][0], "+", '+'*len(master_dict[MI]["Read1"][0])])) #+ = quality score of10
            single_count +=1
        if len(master_dict[MI]["Read2"]) == 1: # if read2 only has 1read thn no sscs or dcs redcan exist
            print>>singles, "\n".join(map(str, ["@"+str(MI)+"-"+str(single_count), master_dict[MI]["Read2"][0], "+", '+'*len(master_dict[MI]["Read2"][0])])) #+ = quality score of10
            single_count +=1
        if len(master_dict[MI]["Read1"]) == 1 and len(master_dict[MI]["Read2"]) == 1:
            #doube singletons included within singleton files, if going include this it needsto be done durring consensus sequence gen
            double_singletons += 1

########Reporting information

end = datetime.datetime.now()
elapsed = end - start

print>>log, "Histogram of number of reads per MI."
for key in sorted(Reads_per_MI_Dict.keys()):
    print>>log, key, "\t", Reads_per_MI_Dict[key]

print>>log, "Accuracy info divided by number MI read depth. Rows = Read Depth, Columns = base position"
print>>log, "\t".join(map(str, ["Read Depth"]+list(sorted(base_positions))))
for read_Depth in Base_discrepancies:
    to_print = [read_Depth]
    for position in list(sorted(base_positions)):
        to_print.append(numpy.mean(Base_discrepancies[read_Depth][position]))
#        print Base_discrepancies[read_Depth][position]
    print>>log, "\t".join(map(str, to_print))

print>>log, "Histogram info on number of discrepent bases:"
print>>log, "\t".join(map(str, ["Read Position", "Unresolved Consensus Seq", "total Consesnsus calls"]))
#for entry in sorted(Unresolved_consensus.keys()):
for position in list(sorted(base_positions)):
    print>>log, position, "\t", Unresolved_consensus[position], "\t", (total_concensus_attempts - Unresolved_consensus[position])

print>>log, "total time elapsed:\t%s" %elapsed
print>>log, "Total Reads:\t%s" %(line_count/2) #/2 because forward and reverse reads counted separate
print>>log, "Dual MI Reads:\t%s" %good_reads
print>>log, "Total MIs:\t%s" %len(master_dict)
print>>log, "DCS reads:\t%s" %dcs_count
print>>log, "SSCS count:\t%s" %sscs_count
print>>log, "single reads:\t%s" %single_count
print>>log, "double singletons:\t%s" %double_singletons
