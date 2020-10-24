# You need to create a password and username here to use with this script:
#   https://public-registry.jbei.org/register
# Note, the script is re-entrant, so if you have an error, just set
#  --initial-chunk to X if the last existing file is part_chunk_X.jsoon

import argparse
from Bio.Seq import Seq
import csv
import json
import os
import pathlib
import time

# How many seconds to wait between queries to not overload requests
query_wait_time = 0

parser = argparse.ArgumentParser(description='download all sequences in JBEI-ICE into the current directory')

parser.add_argument(
    '-e', '--email',
    action='store',
    dest='email',
    required=True,
    type=str,
    help="email (JBEI-username)",
)

parser.add_argument(
    '-p','--password',
    action='store',
    dest='password',
    required=True,
    type=str,
    help="JBEI-ICE password",
)

parser.add_argument(
    '-o', '--output',
    action='store',
    dest='output',
    required=False,
    type=str,
    help="output filename",
    default='JBEI-ICE'
)

parser.add_argument(
    '-c', '--initial-chunk',
    action='store',
    dest='initial_chunk',
    required=False,
    type=int,
    help="start at this chunk",
    default=0
)

args = parser.parse_args()


## Create the folders
output_folder = os.path.join(args.output)
entries_folder = os.path.join(output_folder, 'parts')
sequences_folder = os.path.join(output_folder, 'sequences')


pathlib.Path(output_folder).mkdir(exist_ok=True)
pathlib.Path(entries_folder).mkdir(exist_ok=True)
pathlib.Path(sequences_folder).mkdir(exist_ok=True)


access_token_filename = os.path.join(output_folder, 'access.token.json')
access_token_cmd = "curl -X POST -H \"Content-Type: application/json\" -d '{email:\"" + args.email + "\", password:\"" + args.password + "\"}' \"https://public-registry.jbei.org/rest/accesstokens\" > " + access_token_filename
print(access_token_cmd)
os.system(access_token_cmd)
time.sleep(query_wait_time)
session_id = ''
with open(access_token_filename) as json_file:
    data = json.load(json_file)
    session_id = data['sessionId']
    print(data)
    time.sleep(query_wait_time)

nonempty_page = True
num_entries = 0
num_entries_with_sequences = 0
num_features = 0

entry_page_index = args.initial_chunk
entries_per_page = 100 # Don't set higher than 1000 or website will time out!

features_columns = []
entries_columns = []

while nonempty_page:
    all_entries = []
    all_features = []
    entry_offset = entry_page_index * entries_per_page
    print("Downloading entries: " + str(entry_offset) + "-" + str(entry_offset+entries_per_page-1))

    entry_search_filename = os.path.join(entries_folder, "part_chunk_" + str(entry_page_index) + ".json")
    entry_search_cmd = "curl -X GET -H \"Content-Type: application/json\" -H \"X-ICE-Authentication-SessionId: " + session_id + "\" -H \"Cache-Control: no-cache\" \"https://public-registry.jbei.org/rest/search?q=&limit=" + str(entries_per_page) + "&offset=" + str(entry_offset) + "\" > " + entry_search_filename
    print(entry_search_cmd)
    os.system(entry_search_cmd)
    time.sleep(query_wait_time)

    this_num_entries = 0
    with open(entry_search_filename) as json_file:
        data = json.load(json_file)
        for entry in data['results']:
            this_num_entries = this_num_entries + 1
            num_entries = num_entries + 1
            this_entry = entry['entryInfo']
            all_entries.append(this_entry)
            if (len(entries_columns)==0):
                entries_columns = this_entry.keys()

            this_entry['sequence_length'] = 0
            this_entry['sequence'] = ''

            #Also download the sequence`
            if (this_entry['hasSequence'] == True):
                num_entries_with_sequences = num_entries_with_sequences + 1
                sequence_download_filename = os.path.join(sequences_folder, "part_sequence_" + str(this_entry['id']) + ".json")
                sequence_download_cmd = "curl -X GET -H \"Content-Type: application/json\" -H \"X-ICE-Authentication-SessionId: " + session_id + "\" -H \"Cache-Control: no-cache\" \"https://public-registry.jbei.org/rest/parts/" + str(this_entry['id']) + "/sequence\" > " + sequence_download_filename
                print(sequence_download_cmd)
                os.system(sequence_download_cmd)
                time.sleep(query_wait_time)

                with open(sequence_download_filename) as sequence_json_file:
                    sequence_data = json.load(sequence_json_file)
                    # Some sequences have a mix of RNA and DNA... change to DNA
                    part_sequence = Seq(sequence_data['sequence'].lower().replace('u', 't'))
                    this_entry['sequence_length'] = len(part_sequence)
                    this_entry['sequence'] = str(part_sequence)

                    for feat in sequence_data['features']:
                        #print(feat)

                        # Locations can have multiple pieces (probably for crossing origin?)
                        # Also, sometimes there are no locations!! We set start and end to -1 when this is the case
                        seq_pieces = []
                        for loc in feat['locations']:
                            feat_sequence = part_sequence[(loc['genbankStart']-1):(loc['end'])]
                            if (feat['strand'] == -1):
                                feat_sequence = feat_sequence.reverse_complement()
                            seq_pieces.append(str(feat_sequence))

                        if (feat['strand'] == -1):
                            seq_pieces = reversed(seq_pieces)

                        full_feat_seq = "".join(seq_pieces)

                        add_feat = {}
                        add_feat['part_id'] = this_entry['id']
                        add_feat['id'] = feat['id']
                        add_feat['name'] = feat['name']
                        add_feat['type'] = feat['type']
                        add_feat['start'] = -1
                        add_feat['end'] = -1
                        if (len(feat['locations']) > 0):
                            add_feat['start'] = feat['locations'][0]['genbankStart']
                            add_feat['end'] = feat['locations'][-1]['end']
                        add_feat['strand'] = feat['strand']
                        add_feat['sequence'] = full_feat_seq
                        add_feat['multiple_locations'] =  len(feat['locations']) > 1

                        #print(add_feat)
                        all_features.append(add_feat)

                        if (len(features_columns)==0):
                            features_columns = add_feat.keys()

                        num_features = num_features + 1




    if (this_num_entries == 0):
        nonempty_page = False


    # Write main entry info to a CSV
    outputMatchFile = open(os.path.join(entries_folder,"part_chunk_" + str(entry_page_index) + ".csv"), "w")
    matches_writer= csv.DictWriter(outputMatchFile, extrasaction='ignore', fieldnames=entries_columns)
    matches_writer.writeheader()
    for row in all_entries:
        matches_writer.writerow(row)

    # Write main entry info to a CSV
    outputMatchFile = open(os.path.join(sequences_folder,"features_chunk_" + str(entry_page_index) + ".csv"), "w")
    matches_writer= csv.DictWriter(outputMatchFile, extrasaction='ignore', fieldnames=features_columns )
    matches_writer.writeheader()
    for row in all_features:
        matches_writer.writerow(row)

    #print(all_entries)
    entry_page_index = entry_page_index+1

print("Total entries downloaded      " + str(num_entries))
print("Total entries with sequences  " + str(num_entries_with_sequences))
print("Total features                " + str(num_features))

all_entries = []
all_features = []

entry_page_index = 0
while os.path.exists(os.path.join(entries_folder,"part_chunk_" + str(entry_page_index) + ".csv")):
    inputMatchFile = open(os.path.join(entries_folder,"part_chunk_" + str(entry_page_index) + ".csv"))
    matches_reader= csv.DictReader(inputMatchFile)
    for row in matches_reader:
        all_entries.append(row)
    entry_page_index = entry_page_index + 1

# Write main entry info to a CSV
outputMatchFile = open(os.path.join(output_folder,"all_parts.csv"), "w")
matches_writer= csv.DictWriter(outputMatchFile, extrasaction='ignore', fieldnames=sorted(all_entries[0].keys()) )
matches_writer.writeheader()
for row in all_entries:
    matches_writer.writerow(row)


entry_page_index = 0
while os.path.exists(os.path.join(sequences_folder,"features_chunk_" + str(entry_page_index) + ".csv")):
    inputMatchFile = open(os.path.join(sequences_folder,"features_chunk_" + str(entry_page_index) + ".csv"))
    matches_reader= csv.DictReader(inputMatchFile)
    for row in matches_reader:
        all_features.append(row)
    entry_page_index = entry_page_index + 1

# Write main entry info to a CSV
outputMatchFile = open(os.path.join(output_folder,"all_features.csv"), "w")
matches_writer= csv.DictWriter(outputMatchFile, extrasaction='ignore', fieldnames=sorted(all_features[0].keys()) )
matches_writer.writeheader()
for row in all_features:
    matches_writer.writerow(row)





