#!/usr/bin/env python

# Load packages

import pandas as pd
import numpy as np
import csv
from sklearn.metrics.pairwise import pairwise_distances
from metadata_to_layout import metadata_to_layout


# Parse Arguments

def pw_sim(input, out, isol_meta, gene_meta, run_type, sim_metric, isol_filt, gene_filt, threads, output_matrix):
    
    # Specify run type
    # This set-up was used to avoid reading large files into memory twice if calculating both isolate and gene sim together
    
    if run_type not in ("both", "isolates", "genes"):
        print("ERROR: Invalid run type: must be both, isolates or genes")
        exit(1)
        
    # Load Data
    # specifies data types for import; first column string (gene names), all others bool for similarity calc
    
    col_names = pd.read_csv(args.input, nrows=0, sep = '\t').columns
    types_dict = {col_names[0]: str} # assumes first column is gene/allele names
    types_dict.update({col: bool for col in col_names if col not in types_dict}) # all other columns types set to boolean (requires binary 0/1 input)
    
    try:
    	data = pd.read_csv(args.input, sep = '\t', dtype=types_dict)
    except ValueError:
    	print("Error with input file type. Check file is binary matrix")
    except:
    	print("Error when trying to read input file")
    else:
    	data = pd.read_csv(args.input, sep = '\t', dtype=types_dict)
    
    data.columns = data.columns.str.replace('.gff$', '') # removes .gff from end of isolate names (e.g. from PPanGGoLIN ouputs) - comment out if necessary
    
    # Set file names for write out
    
    if args.out:
        file_out_prefix = args.out
        genes_out = f'{file_out_prefix}_genes_pw_sim.layout'
        isols_out = f'{file_out_prefix}_isols_pw_sim.layout'
        isols_matrix_out = f'{file_out_prefix}_isols_matrix.tsv'
        genes_matrix_out = f'{file_out_prefix}_genes_matrix.tsv'
    else:
        isols_out = 'isols_pw_sim.layout'
        genes_out = 'genes_pw_sim.layout'
        isols_matrix_out = 'isols_matrix.tsv'
        genes_matrix_out = 'genes_matrix.tsv'
    
    # For Isolate-Isolate Comparison
    
    if run_type in ("isolates", "both"):
        
        print("\n-----------------------------------------------\n")
        print("Calculating pairwise distance between isolates...\n")
        
        isols = data.iloc[:,1:]
        isols_trans = isols.transpose()
        isols_jac_sim = pd.DataFrame((1 - pairwise_distances(isols_trans.to_numpy(), metric = sim_metric, n_jobs = args.threads)), index=isols.columns, columns=isols.columns)
	
	# round sim values to 5 dp - comment out if you want full values
        isols_jac_sim = isols_jac_sim.round(5)
        if output_matrix:
            print(" - Saving isolate similarity matrix to file\n")
            isols_jac_sim.to_csv(isols_matrix_out, sep="\t")
        
	# convert to pairwise
        
        isols_to_keep = np.triu(np.ones(isols_jac_sim.shape), k=1).astype('bool').reshape(isols_jac_sim.size)
        isols_jac_pw = isols_jac_sim.stack()[isols_to_keep]
        isols_jac_pw.index.rename(['Isolate1', 'Isolate2'], inplace = True)
        isols_jac_pw = isols_jac_pw.to_frame('jacc_sim').reset_index()
        isols_jac_pw = isols_jac_pw[isols_jac_pw.jacc_sim.ge(args.isol_filt)]
        
        print(" - Isolate pairwise distance calculated\n")
    
        print(" - Printing pairwise list to file...\n")
        isols_jac_pw.to_csv(isols_out, sep = '\t', index = False, header = False, encoding='ISO-8859-1')
        print(" - Isolate pairwise distance file printed\n")
        
        # Add in Metadata for isolates
        if args.isol_meta:
            print("Adding Isolate metadata:\n")
            metadata_to_layout(isols_out, isol_meta, selection=0, run_type="append", verbose=1)
        else:
            print("No isolate metadata file provided.\n")

    # For Gene-Gene Comparison
    
    if run_type in ("genes", "both"):
        
        print("\n-----------------------------------------------\n")
        print("Calculating pairwise distance between genes...\n")
        genes = data
        genes.index = genes[col_names[0]]
        genes = genes.drop([col_names[0]], axis = 1)

        genes_jac_dist = pairwise_distances(genes.to_numpy(), metric = sim_metric, n_jobs = args.threads)
        genes_jac_sim = pd.DataFrame((1 - genes_jac_dist), index=genes.index, columns=genes.index)

	# round sim values to 5 dp - comment out if you want full values
        genes_jac_sim = genes_jac_sim.round(5)

        if output_matrix:
            print(" - Saving gene similarity matrix to file\n")
            genes_jac_sim.to_csv(genes_matrix_out, sep="\t")
	
        # convert to pairwise
        
        genes_to_keep = np.triu(np.ones(genes_jac_sim.shape), k=1).astype('bool').reshape(genes_jac_sim.size)
        genes_jac_pw = genes_jac_sim.stack()[genes_to_keep]
        genes_jac_pw.index.rename(['Gene1', 'Gene2'], inplace = True)
        genes_jac_pw = genes_jac_pw.to_frame('jacc_sim').reset_index()
        genes_jac_pw = genes_jac_pw[genes_jac_pw.jacc_sim.ge(args.gene_filt)]
        
        print(" - Gene pairwise distance calculated\n")

        print(" - Printing pairwise list to file...\n")
        
        genes_jac_pw.to_csv(genes_out, sep = '\t', index = False, header = False, encoding='ISO-8859-1')
        
        # Add in Metadata for genes
        if args.gene_meta:
            print("Adding gene metadata:\n")
            metadata_to_layout(genes_out, gene_meta, selection=0, run_type="append", verbose=1)
        else:
            print("No gene metadata provided\n")
    
    print("-----------------------------------------------\n")
    print("Script finished succesfully!\n")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, required = True, help = "binary tab delimited (.tsv, .Rtab) file of gene presc/absc across pangenome")
    parser.add_argument("-o", "--out", type = str, required = False, default = "", help = "optional prefix for output files")
    parser.add_argument("-m", "--isol_meta", type = str, required = False, help = ".csv isolate metadata")
    parser.add_argument("-g", "--gene_meta", type = str, required = False, help = ".tsv gene metadata, e.g. PIRATE all_alleles.tsv file")
    parser.add_argument("-r", "--run_type", type = str, required = False, default = "both", help = "calculate similarity of isolates (\"isolates\"), genes (\"genes\" or both (\"both\"). Defaults to both")
    parser.add_argument("-s", "--sim_metric", type = str, required = False, default = "jaccard", help = "set metric for similarity calculation. Default: jaccard")
    parser.add_argument("-f", "--isol_filt", type = float, required = False, default = 0.5, help = "optional filter for isolate pairwise similarity")
    parser.add_argument("-e", "--gene_filt", type = float, required = False, default = 0.5, help = "optional filter for gene pairwise similarity")
    parser.add_argument("-t", "--threads", type = int, required = False, default = 1, help = "number of threads to be used")
    parser.add_argument("-x", "--output_matrix", action = 'store_true', required = False, help = "output .tsv similarity matrix as well as list. Will not be filtered")
    args = parser.parse_args()
    
    pw_sim(args.input, args.out, args.isol_meta, args.gene_meta, args.run_type, args.sim_metric, args.isol_filt, args.gene_filt, args.threads, args.output_matrix)

