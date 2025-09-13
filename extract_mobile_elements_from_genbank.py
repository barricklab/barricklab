#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Extract mobile elements from a GenBank file into a multi-FASTA file.")
parser.add_argument("-i", "--input", required=True, help="Input GenBank file")
parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
parser.add_argument("-e", "--element", required=True, help="Only output this element")

args = parser.parse_args()

genbank_file = args.input
fasta_file = args.output
specified_element =  args.element

# Store all extracted mobile element sequences
mobile_elements = []

# Parse the GenBank file
for record in SeqIO.parse(genbank_file, "genbank"):
    for feature in record.features:
        # Check for mobile element tags
        if feature.type == "mobile_element" or "mobile_element_type" in feature.qualifiers:
            # Extract sequence based on feature location
            subseq = feature.extract(record.seq)
            
            if (specified_element):
                if "name" in feature.qualifiers:
                    element_type = feature.qualifiers["name"][0]
                    if (element_type != specified_element):
                        continue;
                else:
                    continue;

            # Build a header
            element_id = record.id + "_" + str(feature.location.start+1) + "_" + str(feature.location.end+1) + "_" + str(feature.location.strand)
            if "name" in feature.qualifiers:
                element_type = feature.qualifiers["name"][0]
                element_id += "_" + element_type
            
            # Create a SeqRecord for this mobile element
            mobile_elements.append(
                SeqRecord(subseq, id=element_id, description="")
            )

# Write all mobile elements to FASTA
SeqIO.write(mobile_elements, fasta_file, "fasta")

print(f"Multi-FASTA with mobile elements written to {fasta_file}")
