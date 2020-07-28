
import csv
from sys import stdout

def metadata_to_layout(layout, metadata, selection, run_type, verbose):
    
    # Check and set run type
    
    if run_type in ("copy"):
        # script to create copy of .layout file path
        # layout_file_out = 
        print("Option currently unavailable. Please create manual copy")
        #print("Copy of layout file will be used")
    
    elif run_type in ("append"):
        print("\nMetadata will be appended to specified .layout file\n")
    
    else: 
        print("ERROR: Invalid run type")
        exit(1)
    
    # Select metadata columns
    
    if selection:
        with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')   # note use of tab-delimited here
            headers = next(csv_reader, None)    # first row of csv assumed to be headers
            with open(selection) as f:
                meta_cols = f.read().splitlines()
                meta_cols_indx = []
                for i in range(0, len(meta_cols)):
                    meta_cols_indx.append(headers.index(meta_cols[i]))
    
    else:
        with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')   # note use of tab-delimited here
            headers = next(csv_reader, None)
            meta_cols_indx = list(range(1, len(headers)))
        if verbose:
            print(" - Note: No metadata columns specified. All will be added")
    
        
    # Count rows   
    
    with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        row_count = sum(1 for row in csv_reader)
    
    if verbose:
        print(" -", row_count, "rows of metadata to be processed")
    
    
    # Append (selected) metadata to layout file
    
    with open(layout, 'a') as out_file:
        with open(metadata, 'r', encoding = "ISO-8859-1") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                for i in meta_cols_indx:     # start at 1 to avoid using 0 index, assuming this is the isolate/gene name 
                    if not row[i] in (""):   # avoids printing blanks
                        print(f'//NODECLASS\t\"{row[0]}\"\t\"{row[i]}\"\t\"{headers[i]}\"', file = out_file)     # row[0] prints isolate/gene name
#                        if verbose:
#                            stdout.write(" - \r%d\%" % (csv_reader.line_num*100/row_count)) ####TODO: Implement proper progress reporter
    
    if verbose:
        print("\nMetadata added succesfully\n")

if __name__ == "__main__":
    import argparse
    
    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--layout", type = str, required = True, help = ".layout file to copy and/or append metadata to")
    parser.add_argument("-m", "--metadata", type = str, required = True, help = ".tsv file of metadata to add. First column must match node names")
    parser.add_argument("-s", "--selection", type = str, required = False, help = ".txt file of list of metadata columns to add. If not provided, all metadata columns will be added")
    parser.add_argument("-r", "--run_type", type = str, required = False, default = "append", help = "copy or append layout file. Default: \"append\"")
    parser.add_argument("-v", "--verbose", action = "store_true", help = "default: on")
    args = parser.parse_args()
    
    metadata_to_layout(args.layout, args.metadata, args.selection, args.run_type, args.verbose)
    


