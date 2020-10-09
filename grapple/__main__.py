import argparse

from subprocess import check_output

#from .py_jaccard_sim import jaccard_sim
#from .py_edges_to_layout import edges_to_layout

def main():
	parser = argparse.ArgumentParser(description='GraPPLE (Graphical Processing for Pangenome Linked Exploration',
					 				prog = 'GraPPLE')

	# required for post_pirate_processing
	parser.add_argument('-t', '--thresholds', required = True, type = str)
	parser.add_argument('-q', '--path_to_PIRATE', required = True, type = str)
	parser.add_argument('-p', '--split_paralogs', action = 'store_true')
	parser.add_argument('-d', '--paralog_dir', default = "with-paralogs")

	# required for jaccard distances
	parser.add_argument('-b', '--binary-data', required = True, type = str)

	# optional for jaccard distances
	parser.add_argument('-m', '--metadata', required = False, type = str)
	parser.add_argument('-g', '--gene-data', required = False, type = str)
	parser.add_argument('-r', '--run-type', required = False, default = "both", type = str)
	parser.add_argument('-n', '--threads', required = False, default = 2, type = int)
	parser.add_argument('-o', '--output', required = False, default = '', type = str)

	# required for edges to layout
	parser.add_argument('-e', '--edges', required = True)
	
	# other
	parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s 0.1.0')

	main = parser.parse_args()

	# 1. post_pirate_processing.sh
	test = check_output(["./grapple/test_script.sh", "test", "output"])
	print(test)

	#check_output(["./post_pirate_processing.sh", main.thresholds, str(main.split_paralogs), main.paralog_dir, main.path_to_PIRATE, str(main.threads)])
	


	print("PIRATE files processed. Now calculating similarity matrices...")

	# 2. jaccard similarity
	# jaccard_sim(input = ,
	# 			output = args.out,
	# 			isol_meta = args.metadata,
	# 			gene_meta = ,
	# 			run_type = args.run_type,
	# 			threads = args.threads)

	# 3. generate_edges.sh

	# 4. edges to layout


	return

if __name__ == "__main__":
	main() 
