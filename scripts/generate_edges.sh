#!\bin\bash

# Author: JDHL
# Basic runner script for adapting PIRATE output for use in GrapPLE/Graphia
# Creates gene presc/absc files using alleles as the gene family category
# Then generates edges files from these
# Heavily relies on PIRATE's excellent adapter scripts, provided in the PIRATE repository
# Please see SionBayliss/PIRATE for more information on these

# error handling
set -oe pipefail

# read variables
while [ "$1" != "" ]; do
        case $1 in
                -i | --input )	    	shift
										input_file=$1
                                        ;;
				-o | --output )			shift
										output_folder=$1
										;;
				-q | --path )			shift
										path=$1
										;;
				-n | --threads )		shift
										threads=$1
										;;
				-h | --help )           printf "\n------------------------------------------------\n\n"
										printf "Ensure you are executing this script within the output directory of your pangenome analysis\n"
										printf -- "-i | input file name - e.g. PIRATE.all_alleles.wp.90.tsv\n"
										printf -- "-o | output folder. Default: pangenome_maps/\n"
										printf -- "-q | path to PIRATE directory - necessary for using PIRATE adapter scripts\n"
										printf -- "-n | number of threads to use. Default: 2\n"
										printf "\n------------------------------------------------\n\n"
										exit 0
                                        ;;
                * )                     printf "Use -h | --help"
                                        exit 1
        esac
        shift
done


# Check Inputs
if test -z "$input_file"; then
	printf "\nERROR: No input file provided"
	exit 1
fi

if test -z "$output_folder"; then
	output_folder="pangenome_maps/"
	if test -d "$output_folder"; then
		printf "\nERROR: Output folder already exists"
		exit 1
	fi
fi
mkdir $output_folder

if ! test -e "$path" || ! test -f $path/scripts/pangenome_graph.pl ; then
    printf "\nERROR: Valid path to PIRATE not provided or is missing necessary pangenome_graph.pl script"
    exit 1
fi

if test -z "$threads"; then
    threads=2
fi

# Set output file names
output_file_tsv=${input_file%.tsv}.allelesAsGeneFamily.tsv

# 1. Generate alleles as gene family files - necessary for edges/layout file on higher thresholds
awk -F '\t' 'BEGIN { OFS = FS } ; { if (NR > 1) { $2 = $1 } } ; {print}' $input_file > $output_folder/$output_file_tsv

# 2. Generate edge files at each threshold specified
$path/scripts/pangenome_graph.pl --input $output_folder/$output_file_tsv --output $output_folder --gffs ./modified_gffs/ --no-cluster --gfa1 

# 3. List files generated (for __main__.py)

declare -a files
for i in ${thr_list//,/ }
do
	files+=$output_folder/$output_file_edges
done

echo "${files[@]}"

