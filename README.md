# Graphical Processing for Pangenome Linked Exploration

#### Support scripts for the visualisation of pangenome analyses in Graphia

These scripts are provided to help users visualise their bacterial pangenome dataset in the network analysis suite [Graphia](https://graphia.app). At present, they should be considered as in beta, so please check your outputs make sense, and raise any issues or bugs.

GraPPLE was initially developed to work with the output of [PIRATE](https://github.com/SionBayliss/PIRATE), though any gene presc/absc matrix from other suitable tools (Roary, Panaroo, PPanGGOLiN) can be used as input for the pairwise similarity script. Currently, synteny graphs from PIRATE and [Panaroo](https://github.com/gtonkinhill/panaroo) are supported (with conversion needed for PIRATE, see below). See also the User Guide for more information.

Contact: j.d.harling-lee [at] roslin.ed.ac.uk

## Dependencies
* Python 3.6
   * [sklearn.metrics](https://scikit-learn.org/stable/install.html)
   * [numpy](https://numpy.org/install/)
   * [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)

## Installation
Currently, each script is run individually, so cloning the repository and running whichever scripts you need is the simplest usage method.
`git clone https://github.com/JDHarlingLee/GraPPLE`

## Usage
#### Pairwise Distances - `pw_similarity.py`
Calculates the pairwise distances between genomes and/or genes from a binary matrix

Example: `python pw_similarity.py -i binary_presc_absc.tsv -o example1 -r "both" -s "jaccard" -f 0.8 -e 0.8 -t 2`

  * Requires the gene presc/absc matrix as a binary, tab separated file (see below for help converting file to binary)
  * `-r` specifies the run_type, and can be set as "isolates", "genes" or "both".
  * Optionally include metadata for genomes or genes (these can also be added later)

#### Add Metadata - `metadata_to_layout.py`
Used to add metadata from a table to a graph in .layout format

Example: `python metadata_to_layout.py -l example1_gene_pw_sim.layout -m gene_info.tsv -s pirate_gene_headers.txt -r "copy"`

  * First column of metadata table must match the names of values in the layout file (the "Node Name")
  * If using the gene information from the pangenome tool (e.g. all_alleles.tsv from PIRATE), remember to use the `-s` variable to specify which columns to add
  * `-s` columns should be a list in .txt format, with one column name per row - e.g. the `pirate_gene_headers.txt` file provided in the main GraPPLE folder
  * Metadata can also be added to a network in Graphia through the GUI, see User Guide

Use `script --help` to see full individual script options.

## PIRATE Specific
These scripts are specific to users of [PIRATE](https://github.com/SionBayliss/PIRATE).

#### Edges to Layout - `edges_to_layout.py` 
Converts PIRATE .edge file to .layout file for load to Graphia

Example: `python edges_to_layout.py -e pangenome.edges -o example_graph` 

  * Calls py_metadata_to_layout.py to add provided metadata to file
  * Note default behaviour is to group genes where they appear twice, as directionality is not currently supported
  * If you are wanting to investigate the synteny graph at a specific id threshold (e.g. 90%), you may need to recreate the graph file, see below.
 
#### Edge File at different threshold
Utilises adapter scripts from PIRATE to recreate the synteny graph at a particular threshold

Example: `python generate_edges.py -i PIRATE.all_alleles.90.tsv -o synteny-graph-90 -p /path/to/PIRATE/`

  * Requires input of a presc/absc matrix at a single threshold (see Other Scripts)
  * Requires path to the installation folder of PIRATE (e.g. pip/conda)
  * This script creates "allelesAsGeneFamilies" files, as it replaces the gene_family field with alleles in the presc/absc matrix

## Other Scripts
Provided in the 'scripts' folder are some other useful scripts, including a general script for subsetting PIRATE output files (post_pirate_processing.sh), and a script to create a simple binary file from pangenome tool ouputs (gene_matrix_to_binary.py). Again, use `script -h` to see options available.

## Acknowledgements
These scripts were initially developed from PIRATE outputs, and we thank Sion Bayliss for his advice and useful discussions.

## Citation
