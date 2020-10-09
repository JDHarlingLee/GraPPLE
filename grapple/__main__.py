import argparse
import os
import subprocess

from .py_jaccard_sim import jaccard_sim
from .py_edges_to_layout import edges_to_layout

def main():
	parser = argparse.ArgumentParser(description='GraPPLE (Graphical Processing for Pangenome Linked Exploration',
					 				prog = 'GraPPLE')

	# required for post_pirate_processing
	parser.add_argument('-t', '--thresholds', required = True, type = str)
	parser.add_argument('-q', '--path_to_PIRATE', required = True, type = str)
	parser.add_argument('-p', '--split_paralogs', action = 'store_true')
	parser.add_argument('-d', '--paralog_dir', default = "with-paralogs")

	# optional for jaccard distances
	parser.add_argument('-m', '--metadata', required = False, type = str)
	parser.add_argument('-g', '--gene-data', required = False, type = str)
	parser.add_argument('-r', '--run-type', required = False, default = "both", type = str)
	parser.add_argument('-n', '--threads', required = False, default = 2, type = int)
	parser.add_argument('-o', '--output', required = False, default = '', type = str)

	# required for edges to layout
	parser.add_argument('-u', '--group-genes', action = 'store_false')
	
	# other
	parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s 0.1.0')

	main = parser.parse_args()

	# 1. post_pirate_processing.sh
	post_pirate_processing = subprocess.run(["bash", "grapple/post_pirate_processing.sh", 
							"-t", main.thresholds,
							"-p", str(main.split_paralogs), 
							"-d", main.paralog_dir, 
							"-q", main.path_to_PIRATE],
							cwd=os.path.dirname(os.path.realpath(__file__)), 
							stdout=subprocess.PIPE, 
							shell=True)

	proc_files = post_pirate_processing.stdout.decode('utf-8')
	proc_files = proc_files.strip("\n").split(" ")

	gene_meta = proc_files[::2]
	gene_bin = proc_files[1::2]

	print("\nPIRATE files processed. Now calculating similarity matrices...\n")

	# 2. jaccard similarity
	for i in range(0, len(gene_bin)):
		jaccard_sim(input = gene_bin[i],
					output = main.out,
					isol_meta = main.metadata,
					gene_meta = gene_meta[i],
					run_type = main.run_type,
					threads = main.threads)

	# 3. generate_edges.sh
	edges_files = subprocess.run(["bash", "grapple/post_pirate_processing.sh", 
							"-t", main.thresholds,
							"-p", str(main.split_paralogs), 
							"-d", main.paralog_dir, 
							"-q", main.path_to_PIRATE], 
							cwd=os.path.dirname(os.path.realpath(__file__)),							
							stdout=subprocess.PIPE, 
							shell=True)

	edge_files = result.stdout.decode('utf-8')
	edge_files = proc_files.strip("\n").split(" ")

	# 4. edges to layout
	#for i in range(0, len(edges_files):
	edges_to_layout(edge_file = edge_files[0],
					file_out = main.output)

	return

if __name__ == "__main__":
	main() 
