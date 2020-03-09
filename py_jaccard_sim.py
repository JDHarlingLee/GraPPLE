"""
@author: jdhl
"""

# Load packages

import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
import csv
import argparse


# Parse Arguments

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type = str, required = True, help = "binary .tsv file of gene presc/absc, e.g. from PIRATE")
parser.add_argument("--out", "-o", type = str, required = False, default = "", help = "optional prefix for output files")
parser.add_argument("--run_type", "-r", type = str, required = False, default = "both", help = "calculate similarity of isolates (\"isolates\"), genes (\"genes\" or both (\"both\"). Defaults to both")
parser.add_argument("--isol_meta", "-m", type = str, required = False, help = ".csv isolate metadata")
parser.add_argument("--gene_meta", "-g", type = str, required = False, help = ".tsv gene metadata, e.g. PIRATE all_alleles.tsv file")
parser.add_argument("--isol_filt", "-f", type = float, required = False, default = 0, help = "optional filter for isolate pairwise similarity")
parser.add_argument("--gene_filt", "-e", type = float, required = False, default = 0, help = "optional filter for gene pairwise similarity")
parser.add_argument("--threads", "-t", type = int, required = False, default = 1, help = "number of threads to be used")
args = parser.parse_args()


# Specify run type

if args.run_type in ("both", "isolates", "genes"):
    run_type = args.run_type
else:
    print("Invalid run_type")
    exit(1)


# Load Data
# specifies data types for import; first column string (gene names), all others bool for jaccard calc

col_names = pd.read_csv(args.input, nrows=0, sep = '\t').columns
types_dict = {'Gene': str}
types_dict.update({col: bool for col in col_names if col not in types_dict})
data = pd.read_csv(args.input, sep = '\t', dtype=types_dict)


# Set file names for write out

file_out_prefix = args.out
genes_out = f'{file_out_prefix}_genes_jacc_dist_pw.txt'
isols_out = f'{file_out_prefix}_isols_jacc_dist_pw.txt'


# For Isolate-Isolate Comparison

if run_type in ("isolates", "both"):
    
    print("------------------------------------------------")
    print("Calculating jaccard distance between isolates...")
    
    isols = data.iloc[:,1:]
    isols_trans = isols.transpose()
    isols_jac_sim = pd.DataFrame((1 - pairwise_distances(isols_trans.to_numpy(), metric = "jaccard", n_jobs = args.threads)), index=isols.columns, columns=isols.columns)
    
    # convert to pairwise
    
    isols_to_keep = np.triu(np.ones(isols_jac_sim.shape), k=1).astype('bool').reshape(isols_jac_sim.size)
    isols_jac_pw = isols_jac_sim.stack()[isols_to_keep]
    isols_jac_pw.index.rename(['Isolate1', 'Isolate2'], inplace = True)
    isols_jac_pw = isols_jac_pw.to_frame('jacc_sim').reset_index()
    isols_jac_pw = isols_jac_pw[isols_jac_pw.jacc_sim.ge(args.isol_filt)]
    
    print("------------------------------------------------")
    print("Isolate pairwise distance calculated")
    
    print("------------------------------------------------")
    print("Printing pairwise list to file...")
    
    isols_jac_pw.to_csv(isols_out, sep = '\t', index = False, header = False, encoding='ISO-8859-1')
    
    # Add in Metadata for isolates
    
    if args.isol_meta:
        print("------------------------------------------------")
        print("Adding isolate metadata....")
        
        isol_meta_in = args.isol_meta
        with open(isols_out, 'a', encoding='ISO-8859-1') as out_file: # a+ appends data to file
            with open(isol_meta_in, 'r', encoding='ISO-8859-1') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                headers = next(csv_reader, None)    # first row of csv assumed to be headers - i.e. metadata categories (e.g. host, CC)
                length = len(headers)
                for row in csv_reader:
                    for i in range(1, length):      # start at 1 to avoid using 0 index, assuming this is the isolate name
                        print(f'//NODECLASS\t\"{row[0]}\"\t\"{row[i]}\"\t\"{headers[i]}\"', file = out_file, )     # row[0] prints isolate name
        print("")
        print("...Isolate metadata added")

    
# For Gene-Gene Comparison

if run_type in ("genes", "both"):
    
    print("---------------------------------------------")
    print("Calculating jaccard distance between genes...")
    genes = data.transpose()
    genes.columns = genes.iloc[0]
    genes = genes.drop(['Gene'])
    genes_trans = genes.transpose() # why are we transposing this twice?!
    genes_jac_dist = pairwise_distances(genes_trans.to_numpy(), metric = "jaccard", n_jobs = args.threads)
    genes_jac_sim = pd.DataFrame((1 - genes_jac_dist), index=genes.columns, columns=genes.columns)
    
    # convert to pairwise
    
    genes_to_keep = np.triu(np.ones(genes_jac_sim.shape), k=1).astype('bool').reshape(genes_jac_sim.size)
    genes_jac_pw = genes_jac_sim.stack()[genes_to_keep]
    genes_jac_pw.index.rename(['Gene1', 'Gene2'], inplace = True)
    genes_jac_pw = genes_jac_pw.to_frame('jacc_sim').reset_index()
    genes_jac_pw = genes_jac_pw[genes_jac_pw.jacc_sim.ge(args.gene_filt)]
    
    print("---------------------------------------------")
    print("Gene pairwise distance calculated")
    
    print("---------------------------------------------")
    print("Printing pairwise list to file...")

    genes_jac_pw.to_csv(genes_out, sep = '\t', index = False, header = False, encoding='ISO-8859-1')

    # Add in Metadata for genes
    
    if args.gene_meta:
        
        print("---------------------------------------------")
        print("Adding gene metadata...")
    
        gene_meta_in = args.gene_meta
        with open(genes_out, 'a', encoding='ISO-8859-1') as out_file: # a+ appends data to file
            with open(gene_meta_in, 'r', encoding='ISO-8859-1') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                headers = next(csv_reader, None)    # first row of csv assumed to be headers - i.e. metadata categories (e.g. host, CC)
                for row in csv_reader:
                    for i in range(1, 20):      # start at 1 to avoid using 0 index, assuming this is the isolate name
                        print(f'//NODECLASS\t\"{row[0]}\"\t\"{row[i]}\"\t\"{headers[i]}\"', file = out_file)     # row[0] prints gene name
                        
        print("")
        print("...Gene metadata added")


print("----------------------------")
print("Script finished succesfully!")



