#!/usr/bin/env python

# simple script to convert a generic gene presc/absc matrix to binary format
# removes metadata columns

import csv
import sys

def file_to_binary(input, output, start_col, delimiter):

    start_col = int(args.start_col)
    input_file = args.input
    out_file = args.output
    delimiter = args.delimiter

    with open(out_file, 'w') as output:
        writer=csv.writer(output, delimiter = '\t', lineterminator = '\n')
        
        with open(input_file, 'r') as i:
        
            tsv_file = csv.reader(i, delimiter = delimiter)
            
            # read genomes from header
            header = next(tsv_file)
            
            # rough check if delimiter is correct
            if len(header)<=1:
                sys.exit("Input has too few columns, check file and delimiter")
            
            # take genomes, add first column (gene/allele name)
            binary_header = header[start_col:]
            binary_header.insert(0, header[0])
            
            # write genome names as first row to output file
            writer.writerow(binary_header)
            
            for line in tsv_file:
                # checks for value (making use of pythons blank = False)
                binary_line = [1 if x else 0 for x in line[start_col:]]
                
                # add first col (gene/allele name) at start of line list
                binary_line.insert(0, line[0])
                
                # check to make sure line is the same length as the header
                if not len(binary_line) == len(binary_header):
                    sys.exit("Line %s has incorrect number of values" % tsv_file.line_num)
                
                # write line to file
                writer.writerow(binary_line)
            

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, required = True, help = "gene presc/absc file, e.g. from PIRATE, Roary")
    parser.add_argument("-o", "--output", type = str, required = True, default = "gene_pa_binary.tsv", help = "name for output file")
    parser.add_argument("--start_col", type = str, required = True, help = "start col of individual genome info - e.g. 20 for PIRATE, 15 for Roary")
    parser.add_argument("--delimiter", type = str, required = False, default = '\t', help = "set input file delimiter (assumes tab by default)")
    args = parser.parse_args()

    file_to_binary(args.input, args.output, args.start_col, args.delimiter)
