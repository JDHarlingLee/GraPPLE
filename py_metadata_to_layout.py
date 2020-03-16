# -*- coding: utf-8 -*-
"""
@author: JDHL
"""

import csv
from sys import stdout
import argparse


# Parse Arguments

parser = argparse.ArgumentParser()
parser.add_argument("--layout" , "-l", type = str, required = True, help = ".layout file to copy and/or append metadata to")
parser.add_argument("--metadata", "-m", type = str, required = True, help = ".tsv file of metadata to add. First column must match node names")
parser.add_argument("--selection", "-s", type = str, required = False, help = ".txt file of list of metadata columns to add. If not provided, all metadata columns will be added")
parser.add_argument("--run_type", "-r", type = str, required = False, default = "copy", help = "copy or append layout file. Default: \"copy\"")
args = parser.parse_args()

layout_file_out = args.layout
metadata_file_in = args.metadata
meta_select_file_in = args.selection


# Check and set run type

if args.run_type in ("copy"):
    # script to create copy of .layout file path
#    layout_file_out = 
    print("Copy of layout file will be used")

elif args.run_type in ("append"):
    print("Metadata will be appended to specified .layout file")

else: 
    print("Invalid run type")
    exit(1)


# Select metadata columns

if args.selection:
    with open(metadata_file_in, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')   # note use of tab-delimited here
        headers = next(csv_reader, None)    # first row of csv assumed to be headers
        length = len(headers)
        with open(meta_select_file_in) as f:
            meta_cols = f.read().splitlines()
            meta_cols_indx = []
            for i in range(0, len(meta_cols)):
                meta_cols_indx.append(headers.index(meta_cols[i]))

else:
    with open(metadata_file_in, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')   # note use of tab-delimited here
        headers = next(csv_reader, None)
        meta_cols_indx = list(range(1, len(headers)))
    print("No metadata columns specified. All will be added")

    
# Count rows   

with open(metadata_file_in, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    row_count = sum(1 for row in csv_reader)

print(row_count, "rows of metadata to be processed")


# Append (selected) metadata to layout file

with open(layout_file_out, 'a') as out_file:
    with open(metadata_file_in, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            for i in meta_cols_indx:     # start at 1 to avoid using 0 index, assuming this is the isolate/gene name 
                if not row[i] in (""):   # avoids printing blanks
                    print(f'//NODECLASS\t\"{row[0]}\"\t\"{row[i]}\"\t\"{headers[i]}\"', file = out_file)     # row[0] prints isolate/gene name
                    stdout.write("\r%d" % (csv_reader.line_num*100/row_count))


print("Script complete")
exit(0)


