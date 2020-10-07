#!\bin\bash

# Author: JDHL
# Basic runner script for adapting PIRATE output for use in GrapPLE/Graphia
# Creates gene presc/absc files using alleles as the gene family category
# Then generates edges files from these
# Heavily relies on PIRATE's excellent adapter scripts, provided in the PIRATE repository
# Please see SionBayliss/PIRATE for more information on these

while [ "$1" != "" ]; do
        case $1 in
                -t | --thr_list )	    shift
										thr_list=$1
                                        ;;
                -p | --paralogs )		shift
										paralogs=1
										;;
				-d | --paralog_dir )	shift
										paralog_dir=$1
										;;
				-q | --path )			shift
										path=$1
										;;
				-n | --threads )		shift
										threads=$1
										;;
				-h | --help )           echo "------------------------------------------------\n"
										echo "Ensure you are executing this script within the output directory of your pangenome analysis\n"
                                        echo "Runs from files PIRATE.all_alleles.wp.{thr}.tsv - at each threshold specified"
										echo "-t | threshold list"
										echo "-p | include paralogs or not. Default: off"
										echo "-d | paralog_directory name - new files are created to avoid overwriting originals"
										echo "-q | path to PIRATE directory - necessary for using PIRATE adapter scripts"
										echo "-n | number of threads to use\n"
										echo "------------------------------------------------"
										exit 0
                                        ;;
                * )                     echo "Use -h | --help"
                                        exit 1
        esac
        shift
done


# Check Inputs
if test -z "$thr_list"; then
        echo "ERROR: No thresholds supplied"
        exit 1
fi

if test -z "$paralogs"; then
		paralogs=0
fi

if test -z "$path"; then
        echo "ERROR: No path to PIRATE provided"
        exit 1
fi

if test -z "$threads"; then
        threads=2
fi

for i in ${thr_list//,/ }
do
    if test -f PIRATE.all_alleles.wp.${i}.tsv; then
        echo "PIRATE.all_alleles.wp.${i}.tsv not found"
        exit 1
    fi
done

# 1. Generate alleles as gene family files - necessary for edges/layout file on higher thresholds
for i in ${thr_list//,/ }
do
	awk -F '\t' 'BEGIN { OFS = FS } ; { if (NR > 1) { $2 = $1 } } ; {print}' PIRATE.all_alleles.wp.${i}.tsv > pangenome_maps/PIRATE.all_alleles.wp.${i}.allelesAsGeneFamily.tsv 
done


# 2. Generate edge files at each threshold specified
for i in ${thr_list//,/ }
do 
   $path/scripts/pangenome_graph.pl --input with-paralogs/PIRATE.all_alleles.wp.${i}.allelesAsGeneFamily.tsv --output pangenome_maps/thr_${i}/ --gffs modified_gffs --no-cluster --gfa1 
done