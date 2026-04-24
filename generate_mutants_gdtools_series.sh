#!/usr/bin/env bash
#Author: Ira Zibbu
#Last update: 2024-02-19
#Description: A script that accepts a list of GenomeDiff files and a reference and sequentially applied the GenomeDiffs
#Usage generate_mutants_gdtools_series.sh <reference.fasta> [file1.gd] [file2.gd] ...
# A conda env with breseq should be active

# Check if at least one filename is provided
if [ "$#" -eq 0 ]; 
then
    echo "Usage: $0 <reference.fasta> [file1.gd] [file2.gd] ..."
    exit 1
fi

fasta_reference=$1
reference_prefix=$(echo ${fasta_reference} | sed "s;.fasta;;")
echo "Reference is ${fasta_reference}"
shift # shift the index of positional args. Now the the first gd file is in $1
# Loop through each provided file 

counter=0 # counter for the names of the output files
for file in "$@"; do

if [ $counter -eq 0 ];
then
gdtools APPLY -o "${reference_prefix}.1.fasta" -f fasta -r ${fasta_reference} $file
((counter++))
continue 1
fi

input_filename="${reference_prefix}.${counter}.fasta"
echo ${input_filename}
((counter++))
output_filename="${reference_prefix}.${counter}.fasta"

gdtools APPLY -o ${output_filename} -f fasta -r ${input_filename} $file
done
