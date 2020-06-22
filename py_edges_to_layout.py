
import csv

def edges_to_layout(edges, file_out, filter, verbose):
    written = 0
    removed = 0
    
    if verbose:
        print("\n--------------------------------------\n")
        print("Converting edge file to .layout file:\n")
    
    if file_out:
        file_out = file_out + ".layout"
    else:
        file_out = "pangenome.layout"
    
    with open(file_out, 'w') as outfile:
        
        with open(edges, 'r') as edges:
            
            edge_reader = csv.reader(edges, delimiter = '\t')
            
            for row in edge_reader:
                
                if int(row[2]) > filter:
                    print(f'"{row[0]}\"\t\"{row[1]}\"\t\"{row[2]}\"', file = outfile)
                    written += 1
                
                else:
                    removed += 1
                                              
            if verbose:
                print(" -", written, "edges written, with", removed, "edges removed\n")
                print("Edge file succesfully converted to", file_out, "\n")
                print("-------------------------------\n")
            
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--edges", type = str, required = True, help = "pangenome.edges file from PIRATE")
    parser.add_argument("-o", "--file_out", type = str, required = False, default = "", help = "output prefix")
    parser.add_argument("-f", "--filter", type = int, required = False, default = 0, help = "filter edges below this value. Must be integer value. Default: 0")
    parser.add_argument("-v", "--verbose", action = "store_false", help = "default: on")
    args = parser.parse_args()
    
    edges_to_layout(args.edges, args.file_out, args.filter, args.verbose)



