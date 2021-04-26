#!/usr/bin/env python


import csv
import os
from sys import stdout
from shutil import copyfile

def metadata_to_layout(layout, metadata, selection, run_type, verbose):
    
    if not os.path.isfile(layout):
        print("\nERROR: .layout file not found.")
        exit(1)

    # Check and set run type
    if run_type in ("c", "cp", "copy"):
        file_base = os.path.splitext(layout)[0]
        file_copy_name = str(file_base + '-wMeta.layout')
        copyfile(layout, file_copy_name)
        layout = file_copy_name
        if verbose:
            print("\nMetadata will be added to the copied file %s\n" % layout)

    elif run_type in ("a", "ap", "append"):
        if verbose:
            print("\nMetadata will be appended directly to %s\n" % layout)
    
    else: 
        print("\nERROR: Invalid run type\n")
        exit(1)

    # Check and set metadata delimiter
    if metadata.endswith('.csv'):
        metadata_delim = ','
    elif metadata.endswith('.tsv'):
        metadata_delim = '\t'
    elif metadata.endswith('.Rtab'):
        metadata_delim = '\t'
    else:
        print("ERROR: Invalid metadata file type. Please provide in either .csv or .tsv format")
        exit(1)
    
    # Select metadata columns
    if selection:
        with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=metadata_delim)   # note use of tab-delimited here
            headers = next(csv_reader, None)    # first row of csv assumed to be headers
            with open(selection) as f:
                meta_cols = f.read().splitlines()
                meta_cols_indx = []
                for i in range(0, len(meta_cols)):
                    meta_cols_indx.append(headers.index(meta_cols[i]))
                if verbose:
                    print(" - %s metadata columns will be added" % len(meta_cols))
    
    else:
        with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=metadata_delim)   # note use of tab-delimited here
            headers = next(csv_reader, None)
            meta_cols_indx = list(range(1, len(headers)))
            if verbose:
                print(" - No columns specified. All %s will be added" % len(headers))
    
        
    # Count rows   
    
    with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=metadata_delim)
        row_count = sum(1 for row in csv_reader)
        if verbose:
            print(" - %s rows of metadata will be added" % row_count)
    
    
    # Append (selected) metadata to layout file
    
    with open(layout, 'a') as out_file:
        out_file.write('\n') # ensures metadata is added beginning on new line
        with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=metadata_delim)
            for row in csv_reader:
                for i in meta_cols_indx:     # start at 1 to avoid using 0 index, assuming this is the isolate/gene name 
                    if not row[i] in (""):   # avoids printing blanks
                        print(f'//NODECLASS\t\"{row[0]}\"\t\"{row[i]}\"\t\"{headers[i]}\"', file = out_file)     # row[0] prints isolate/gene name

    if verbose:
        print("\nAll Metadata added succesfully\n")

if __name__ == "__main__":
    import argparse
    
    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--layout", type = str, required = True, help = ".layout file to copy and/or append metadata to")
    parser.add_argument("-m", "--metadata", type = str, required = True, help = ".tsv file of metadata to add. First column must match node names")
    parser.add_argument("-s", "--selection", type = str, required = False, help = ".txt file of list of metadata columns to add. If not provided, all metadata columns will be added")
    parser.add_argument("-r", "--run_type", type = str, required = False, default = "copy", help = "copy or append layout file. Default: \"copy\"")
    parser.add_argument("-v", "--verbose", action = "store_false", help = "default: on")
    args = parser.parse_args()
    
    metadata_to_layout(args.layout, args.metadata, args.selection, args.run_type, args.verbose)
    


