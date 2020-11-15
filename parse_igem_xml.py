#This script can parse the XML download of the iGEM Registry database available here:
# http://parts.igem.org/partsdb/download.cgi?type=parts
# 20 November 2019 version
#
# After downloading this, you need to remove "gremlin" characters from it, most of which occur in the encoded sha1 key entries
# I did this by using TextWrangler's "Zap Gremlins" command and deletig those characters
#
# Save that input file as: xml_parts_gremlins_zapped.xml
# Run this script from that directory
# $ python3 parse_igem_xml.py
# The output is in "xml_parts.csv"

import xmltodict
from Bio.Seq import Seq
import csv

with open('xml_parts_gremlins_zapped.xml', encoding = "utf8") as fd:
    doc = xmltodict.parse(fd.read())

parts = doc["mysqldump"]["database"]["table_data"][0]
parts_seq_features = doc["mysqldump"]["database"]["table_data"][1]

# Step 1: Create a new dictionary that allows looking up part sequencew and info by their part_id

# Entries in the parts database look like this:
# 	<row>
# 		<field name="part_id">19713</field>
# 		<field name="ok">1</field>
# 		<field name="part_name">BBa_K346009</field>
# 		<field name="short_desc">Constutitive promoter(BBa_J23103)+RBS+MerR(BBa_K346001)</field>
# 		<field name="description">Constutitive promoter BBa_J23103+MerR
# </field>
# 		<field name="part_type">Generator</field>
# 		<field name="author">Qianzhu Wu</field>
# 		<field name="owning_group_id">9</field>
# 		<field name="status">Planning</field>
# 		<field name="dominant">1</field>
# 		<field name="informational">0</field>
# 		<field name="discontinued">0</field>
# 		<field name="part_status"></field>
# 		<field name="sample_status">It's complicated</field>
# 		<field name="p_status_cache"></field>
# 		<field name="s_status_cache"></field>
# 		<field name="creation_date">2010-10-18</field>
# 		<field name="m_datetime">2018-11-06 15:23:13</field>
# 		<field name="m_user_id">0</field>
# 		<field name="uses">0</field>
# 		<field name="doc_size">2122</field>
# 		<field name="works">Works</field>
# 		<field name="favorite">0</field>
# 		<field name="specified_u_list">_7145_19320_</field>
# 		<field name="deep_u_list">_7145_19320_</field>
# 		<field name="deep_count">2</field>
# 		<field name="ps_string" xsi:nil="true" />
# 		<field name="scars"></field>
# 		<field name="default_scars"></field>
# 		<field name="owner_id">5916</field>
# 		<field name="group_u_list">_457_</field>
# 		<field name="has_barcode">0</field>
# 		<field name="notes">Constutitive promoter BBa_J23103+MerR
# </field>
# 		<field name="source">Constutitive promoter BBa_J23103+MerR
# </field>
# 		<field name="nickname"></field>
# 		<field name="categories"></field>
# 		<field name="sequence">ctgatagctagctcagtcctagggattatgctagctactagagaaagaggagaaatactagatggaaaataatttggaaaacctgaccattggcgtttttgccaaggcggccggggtcaacgtggagacaatccgcttctatcagcgcaagggcctgttgcgggaaccggacaagccttacggcagcatccgccgctatggggaggcggacgtggttcgggtgaaattcgtgaaatcggcacagcggctggggttcagtctggacgagattgccgagctgttgcggctcgacgatggcacccactgcgaggaggccagcagcctggccgaacacaagctcaaggacgtgcgcgagaagatggccgacttggcgcgcatggaaaccgtgctgtctgaactcgtgtgcgcctgccatgcacgaaaggggaatgtttcctgcccgttgatcgcgtcactacagggcgaagcaggcctggcaaggtcagctatgccttag</field>

# 		<field name="sequence_update">10</field>
# 		<field name="seq_edit_cache">&lt;script src='http://parts.igem.org/cgi/partsdb/seq_edit/all.js'&gt;&lt;/script&gt; &lt;DIV id='sequencePaneDiv' style='clear:both'&gt; &lt;INPUT type='hidden' id='new_dna_format' name='new_dna_format' value='' /&gt; &lt;INPUT type='hidden' id='selection_start' name='selection_start' value='14' /&gt; &lt;INPUT type='hidden' id='selection_end'   name='selection_end'   value='0' /&gt;&lt;/DIV&gt; &lt;script&gt; var sequence = new String('ctgatagctagctcagtcctagggattatgctagctactagagaaagaggagaaatactagatggaaaataatttggaaaacctgaccattggcgtttttgccaaggcggccggggtcaacgtggagacaatccgcttctatcagcgcaagggcctgttgcgggaaccggacaagccttacggcagcatccgccgctatggggaggcggacgtggttcgggtgaaattcgtgaaatcggcacagcggctggggttcagtctggacgagattgccgagctgttgcggctcgacgatggcacccactgcgaggaggccagcagcctggccgaacacaagctcaaggacgtgcgcgagaagatggccgacttggcgcgcatggaaaccgtgctgtctgaactcgtgtgcgcctgccatgcacgaaaggggaatgtttcctgcccgttgatcgcgtcactacagggcgaagcaggcctggcaaggtcagctatgccttag');        var seqFeatures = new Array( ['brick',44,496,'K346001', 1], ['brick',1,35,'J23103', 1]); var subParts = new Array (  new Part (7145, 'J23103', '', 'http://parts.igem.org/images/partbypart/icon_regulatory.png'),  new Part (19320, 'K346001', '', 'http://parts.igem.org/images/partbypart/icon_translational_unit.png')); var Format = '_subparts_'; var PrimaryPartName = 'BBa_K346009'; var PrimaryPartID = '19713'; var Selection_Start = 0; var Selection_End = 0; showSeqFeatures(false);  &lt;/script&gt;&lt;div style='position:relative;clear:both;width:100%'&gt;&lt;div style=''&gt;
# &lt;STYLE type='text/css'&gt;
# .compatibility_div ul,
# .compatibility_div li {
# display: inline;
# }
# .compatibility_div li {
# position: relative;
# padding-top: 2px;
# padding-left:4px;
# padding-right:3px;
# margin-right:2px;
# margin-bottom: 5px;
# }
# .compatibility_div .box {
# top:   35px;
# width: 200px;
# left:  0px;
# }
# .compatibility_div .box_white {
# border: 1px solid gray;
# background-color: white;
# }
# .compatibility_div .box_red {
# border: 1px solid #dd6666;
# background-color: #ffcccc;
# background-image: url('http://parts.igem.org/images/red_not/red_box.png');
# background-repeat: none;
# }
# .compatibility_div .box_green {
# border: 1px solid #44ee44;
# background-color: #aaffaa;
# }


# &lt;/STYLE&gt;&lt;DIV class='compatibility_div' style='display:inline'&gt;
# Assembly Compatibility:&lt;UL&gt;

# &lt;LI class='boxctrl box_green'&gt;10&lt;DIV class='box'&gt;&lt;DIV&gt;COMPATIBLE WITH RFC[10]&lt;/div&gt;&lt;/div&gt;&lt;/LI&gt;
# &lt;LI class='boxctrl box_red'&gt;12&lt;DIV class='box'&gt;&lt;DIV style='border-bottom:1px solid gray;'&gt;INCOMPATIBLE WITH RFC[12]&lt;/div&gt;Illegal NheI site found at 7&lt;BR&gt;Illegal NheI site found at 30&lt;BR&gt;&lt;/div&gt;&lt;/LI&gt;
# &lt;LI class='boxctrl box_green'&gt;21&lt;DIV class='box'&gt;&lt;DIV&gt;COMPATIBLE WITH RFC[21]&lt;/div&gt;&lt;/div&gt;&lt;/LI&gt;
# &lt;LI class='boxctrl box_green'&gt;23&lt;DIV class='box'&gt;&lt;DIV&gt;COMPATIBLE WITH RFC[23]&lt;/div&gt;&lt;/div&gt;&lt;/LI&gt;
# &lt;LI class='boxctrl box_green'&gt;25&lt;DIV class='box'&gt;&lt;DIV&gt;COMPATIBLE WITH RFC[25]&lt;/div&gt;&lt;/div&gt;&lt;/LI&gt;
# &lt;LI class='boxctrl box_green'&gt;1000&lt;DIV class='box'&gt;&lt;DIV&gt;COMPATIBLE WITH RFC[1000]&lt;/div&gt;&lt;/div&gt;&lt;/LI&gt;
# &lt;/UL&gt;
# &lt;/DIV&gt;&lt;/div&gt;&lt;/div&gt;</field>
# 		<field name="review_result" xsi:nil="true" />
# 		<field name="review_count">0</field>
# 		<field name="review_total">0</field>
# 		<field name="flag" xsi:nil="true" />
# 		<field name="sequence_length">496</field>
# 		<field name="temp_1">0</field>
# 		<field name="temp_2" xsi:nil="true" />
# 		<field name="temp_3" xsi:nil="true" />
# 		<field name="temp4" xsi:nil="true" />
# 		<field name="rating">0</field>
# 	</row>

part_id_dict = {}
num_parts = 0
num_complete_parts = 0

for row in parts["row"]:
    num_parts = num_parts+1
    #print(row)

    # Whatever fields you add to this will be copied from the XML
    new_entry = { 'part_id':'', 'part_name':'', 'sequence':'', 'short_desc':'' }

    for i in range(len(row["field"])):
        field = row["field"][i]

        if (field["@name"] in new_entry.keys()):
            if ("#text" in field):
                new_entry[field["@name"]] = field["#text"]
    
    # Check for incomplete entries that would mess up the dictionary
    if (new_entry['part_id'] == ''):
        continue

    num_complete_parts = num_complete_parts + 1
    part_id_dict[new_entry['part_id']] = new_entry

#print(part_id_dict)

print("Number of parts (complete): " + str(num_parts) + " (" + str(num_complete_parts) + ")")

# Step 2: find the features annotated in those sequences.

# Entries in the 'parts_seq_features' table look like this:
# 	<row>
# 		<field name="feature_id">2015460</field>
#		<field name="feature_type">binding</field>
#		<field name="start_pos">89</field>
#		<field name="end_pos">126</field>
#		<field name="label">CAP binding site</field>
#		<field name="part_id">14520</field>
#		<field name="type" xsi:nil="true" />
#		<field name="label2" xsi:nil="true" />
#		<field name="mark">0</field>
#		<field name="old">0</field>
#		<field name="reverse">0</field> 
#	</row>


output_table = []
num_parts_seq_features = 0
num_complete_parts_seq_features = 0

for row in parts_seq_features["row"]:
    #print(row)

    num_parts_seq_features = num_parts_seq_features + 1

    # Whatever fields you add to this will be copied from the XML
    this_seq_feature = { 'part_id':'', 'feature_id':'', 'feature_type':'', 'start_pos':'', 'end_pos':'', 'label':'', 'reverse':'', 'sequence':'' }

    for i in range(len(row["field"])):
        field = row["field"][i]

        if (field["@name"] in this_seq_feature.keys()):
            if ("#text" in field):
                this_seq_feature[field["@name"]] = field["#text"]
    
    # Now add the sequence if enough information was provided
    # Is the main information provided?
    if (this_seq_feature['label'] == '') or (this_seq_feature['part_id'] == '') or (this_seq_feature['reverse'] == '') or (this_seq_feature['start_pos'] == '') or (this_seq_feature['end_pos'] == ''):
        print("Skipping incomplete seq_feature")
        print(this_seq_feature)
        print(row)
        continue
    
    #Did we find the part_id?
    if ( this_seq_feature['part_id'] not in part_id_dict.keys() ):
        print("Skipping seq_feature due to not finding referenced part_id")
        print(this_seq_feature)
        print(row)
        continue

    this_seq_feature['short_desc'] = part_id_dict[this_seq_feature['part_id']]['short_desc']
    this_seq_feature['part_name'] = part_id_dict[this_seq_feature['part_id']]['part_name']

    #Assign sequence and splice out the right part
    this_seq_feature['sequence'] = part_id_dict[this_seq_feature['part_id']]['sequence']
    this_seq_feature['sequence'] = this_seq_feature['sequence'][(int(this_seq_feature['start_pos'])-1):int(this_seq_feature['end_pos'])]


    num_complete_parts_seq_features = num_complete_parts_seq_features+1

    if (this_seq_feature['reverse'] == "1"):
        #(this_seq_feature)
        temp_seq = Seq(this_seq_feature['sequence'])
        this_seq_feature['sequence'] = str(temp_seq.reverse_complement())
        #print(this_seq_feature)

    output_table.append(this_seq_feature)

#print(output_table)
print("Number of parts (complete): " + str(num_parts) + " (" + str(num_complete_parts) + ")")

outputMatchFile = open("xml_parts.csv", "w")
matches_writer= csv.DictWriter(outputMatchFile, fieldnames=sorted(output_table[1].keys()) )
matches_writer.writeheader()
for row in output_table:
    matches_writer.writerow(row)