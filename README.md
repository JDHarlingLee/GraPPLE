# Graphical Processing for Pangenome Linked Exploration

#### Support scripts for the visualisation of pangenome analyses in Graphia

The below scripts take a bacterial pangenome, calculate similarity matrices between both isolates and genes, and output these for use in the network analysis suite Graphia (https://graphia.app). 
There is also a script for converting the synteny map from PIRATE for use in Graphia (.gfa to .layout).

Currently, the only pangenome tools fully supported is PIRATE (https://github.com/SionBayliss/PIRATE), though in theory any gene presc/absc matrix from another suitable can be used as input.

## Usage
<i> Full Runner Script in progress </i>

1. __post_pirate_processing.sh__ | Convert standard PIRATE output to format necessary for later scripts.
  * Recreates 'all_alleles' files, and can add in paralogs
  * Converts all_alleles to binary files
  * Can be run on all thresholds, or a subset (e.g. "90, 95") to reduce file size

2. __py_jaccard_sim.py__ | Calculates jaccard distance between isolates, genes
  * Requires binary file input
  * Optional include metadata on isolates or genes (can be added later)

3. __py_edges_to_layout.py__ | Converts edges file to layout file for load to Graphia
  * Calls py_metadata_to_layout.py to add provided metadata to file
  * Note default behaviour to group genes as directionality is not currently supported
 
4. __py_metadata_to_layout.py__ | Script to add metadata in either .tsv or .csv format to a .layout file for load to Graphia

Use "--help" to see full individual script options.

## Acknowledgements
These scripts rely heavily on the excellent adapter scripts from PIRATE, and we thank Sion Bayliss also for his advice and useful discussions.

## Citation
