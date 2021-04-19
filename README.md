# Graphical Processing for Pangenome Linked Exploration

#### Support scripts for the visualisation of pangenome analyses in Graphia

GraPPLE allows a user to take their bacterial pangenome dataset, calculate similarity matrices between both isolates and genes, and output these in a format suitable for use in the network analysis suite Graphia (https://graphia.app).

GraPPLE was initially developed to work from the output of PIRATE (https://github.com/SionBayliss/PIRATE), though any gene presc/absc matrix from other suitable tools (Roary, Panaroo, PPanGGOLiN) can be used as input for the pairwise similarity script. See User Guide for more information.

## Usage
1. __py_jaccard_sim.py__ | Calculates jaccard similarity between each pair of genomes and each pair of genes
  * Requires binary, tab separated file input
  * Optionally include metadata for genomes or genes (can be added later)

2. __py_edges_to_layout.py__ | Converts PIRATE .edge file to .layout file for load to Graphia
  * Calls py_metadata_to_layout.py to add provided metadata to file
  * Note default behaviour to group genes as directionality is not currently supported
 
3. __py_metadata_to_layout.py__ | Script to add metadata in either .tsv or .csv format to a .layout file for load to Graphia

Use "--help" to see full individual script options.

## Other Scripts
Provided in the 'scripts' folder are some other useful scripts for compatability between tools, such as a script to create a simple binary file from pangenome tool ouputs.

## Acknowledgements
These scripts were initially developed from PIRATE outputs, and we thank Sion Bayliss for his advice and useful discussions.

## Citation
