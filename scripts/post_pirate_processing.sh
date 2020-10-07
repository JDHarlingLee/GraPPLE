#!\bin\bash

# Author: JDHL
# Basic runner script for adapting PIRATE output for use in GrapPLE/Graphia
# Recreates some files, and splits out files by threshold for calculation of Jaccard Similarity
# Can include paralogs or not (default)
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
										echo "-t | threshold list (those used in PIRATE run, or subset thereof"
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

if [[ $paralogs = 1 ]]; then
	if test -z "$paralog_dir"; then
		paralog_dir="with-paralogs"
	fi
	if test -d "$paralog_dir"; then
		echo "ERROR: paralog output directory already exists. Please choose a different directory name. Exiting to avoid overwriting previous data"
		exit 1
	fi
fi


# For debugging:
#echo $thr_list
#echo $paralogs
#echo $path
#echo $threads
#echo $paralog_dir
#echo "high and low set later in script"
#exit


# 1. Recreate PIRATE.all_alleles.tsv

if test -f ./PIRATE.all_alleles.tsv; then
	echo "PIRATE.all_alleles.tsv already exists"
else 
	$path/scripts/link_clusters_runner.pl -l ./loci_list.tab -t $thr_list -o ./ -c ./co-ords/ -parallel $threads --all-alleles
fi

# 2. Recreate split_paralogs_loci.tab

if test -f ./split_paralog_loci.tab; then
	echo "split_paralogs_loci.tab exists"
else 
	$path/scripts/split_paralogs_runner.pl -p ./loci_paralog_categories.tab -l ./loci_list.tab -o ./ -t $threads
fi

# 3. Create paralog files

if [[ $paralogs = 1 ]]; then
	mkdir ${paralog_dir}/
	${path}/scripts/link_clusters_runner.pl -l ./loci_list.tab -l ./split_paralog_loci.tab -t $thr_list -o ./${paralog_dir}/ -c ./co-ords/ --paralogs ./loci_paralog_categories.tab -e ./paralog_clusters.tab --parallel $threads --all-alleles
	wp="wp."
	mv ${paralog_dir}/PIRATE.all_alleles.tsv ./PIRATE.all_alleles.wp.tsv
else
	wp=""
fi

# 4. Split all_alleles by thresholds

if test -z $high; then
	high=$(cut -f7 ./PIRATE.all_alleles.tsv | sort -n | tail -1)
fi

if test -z $low; then
	low=$(cut -f7 ./PIRATE.all_alleles.tsv | sort -n | head -2 | tail -1)
fi

for i in ${thr_list//,/ }
do
	awk -F"\t" -v thr=$i -v max=$high -v min=$low 'NR==1 {print $0}; $5==thr && $7<max && $7>min {print $0}' PIRATE.all_alleles.${wp}tsv > PIRATE.acc_alleles.${i}.${wp}tsv 
done

# 5. Convert PIRATE files to binary format

for i in ${thr_list//,/ }
do
	$path/tools/convert_format/PIRATE_to_Rtab.pl -i ./PIRATE.acc_alleles.${i}.${wp}tsv -o ./PIRATE.acc_alleles.${i}.${wp}binary.tsv --low 0 --high 1 
done













