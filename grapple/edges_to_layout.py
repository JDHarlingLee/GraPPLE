#!/usr/bin/env python


import pandas as pd 

def edges_to_layout(edge_file, file_out, edge_filter, group_genes, verbose):
    
    if verbose:
        print("\n--------------------------------------\n")
        print("Converting edge file to .layout file:\n")

    if file_out:
        file_out = file_out + ".layout"
    else:
        file_out = "pangenome.layout"
        
# Processing Edge File
    
    with open(file_out, 'w') as outfile:
        
        df = pd.read_csv(edge_file, sep = '\t', header = None)
        
        df.columns = ['gene1', 'gene2', 'edge_weight']
        
        total = len(df)
                
        # remove direction from genes, as "node directionality" is not supported by Graphia
        df = df.replace({'-':''}, regex = True) 
        
        # edges are therefore a sum of connections between genes (e.g. the visual is a simplified representation of gene neighbourhoods)
        if group_genes:
            df = df.groupby(["gene1", "gene2"]).sum()
            grouped = total - len(df)
        
        df = df.reset_index()
        
        df = df[df.edge_weight > edge_filter]
        
        if group_genes:
            filtered = total - grouped - len(df)
        
        df.to_csv(outfile, sep = '\t', index = False, header = False)
        
        written = len(df)
        
        if verbose:
            if group_genes:
                print(" -", grouped, "edges grouped")
                print(" -", filtered, "edges filtered")
            print(" -", written, "edges written\n")
            print("Edge file succesfully converted to", file_out, "\n")
            print("-------------------------------\n")
            

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--edge_file", type = str, required = True, help = "pangenome.edges file from PIRATE")
    parser.add_argument("-o", "--file_out", type = str, required = False, default = '', help = "output prefix")
    parser.add_argument("-f", "--edge_filter", type = int, required = False, default = 0, help = "filter edges below this value. Must be integer value. Default: 0")
    parser.add_argument("-g", "--group_genes", action = "store_false", help = "Default: On. Removes directionality of genes to group together in Graphia. See docs for more information on this behaviour")
    parser.add_argument("-v", "--verbose", action = "store_false", help = "Default: On")
    args = parser.parse_args()
    
    edges_to_layout(args.edge_file, args.file_out, args.edge_filter, args.group_genes, args.verbose)



