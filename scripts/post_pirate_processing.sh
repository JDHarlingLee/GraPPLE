#!/usr/bin/bash

# Author: JDHL
# Basic runner script for adapting the general PIRATE output for use in GraPPLE/Graphia
# Recreates some files, and splits out files by identity threshold
# Can include paralogs or not (default)
# Heavily relies on PIRATE's excellent adapter scripts, provided in the PIRATE repository
# Please see SionBayliss/PIRATE for more information on these

# read variables
while [ "$1" != "" ]; do
        case $1 in
                -t | --thr_list )	shift
					thr_list=$1
                                        ;;
	        -p | --paralogs )	shift
					paralogs=1
					;;
		-d | --grapple_dir )	shift
					grapple_dir=$1
					;;
		-q | --path )		shift
					path=$1
					;;
		-n | --threads )	shift
					threads=$1
					;;
		-h | --help )           printf "\n------------------------------------------------\n"
					printf "\nThis is a script to adapt the standard output from PIRATE for use in GraPPLE/Graphia"
					printf "\nYou may prefer to run some of these actions individually, but this is a wrapper script for convenience\n"
					printf "\nEnsure you are executing this script within the output directory of your pangenome analysis\n"
					printf "\n-t | threshold list (those used in PIRATE run, or a subset thereof)"
					printf "\n-p | include paralogs or not (1 = yes, 0 = no). Default: no"
					printf "\n-d | GraPPLE analysis directory name - a directory where new files can be created to avoid overwriting originals. Default: GraPPLE"
					printf "\n-q | path to PIRATE directory - necessary for using PIRATE adapter scripts"
					printf "\n-n | number of threads to use. Default: 2\n"
					printf "\n------------------------------------------------\n"
					exit 0
                                        ;;
                * )                     printf "\nUse -h | --help"
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

if ! test -e "$path" || ! test -f $path/scripts/link_clusters_runner.pl || ! test -f $path/scripts/split_paralogs_runner.pl || ! test -f $path/tools/convert_format/PIRATE_to_Rtab.pl ; then
        echo "ERROR: Valid path to PIRATE not provided or is missing necessary adapter scripts"
        exit 1
fi

if test -z "$threads"; then
        threads=2
fi

if [[ $paralogs -eq 1 ]]; then
	if test -z "$grapple_dir"; then
		grapple_dir="GraPPLE"
	fi
	if test -d "$grapple_dir"; then
		echo "ERROR: GraPPLE analysis directory (${grapple_dir}) already exists. Please choose a different directory name. Exiting to avoid overwriting previous data"
		exit 1
	fi
fi


#For debugging:
#echo $thr_list
#echo $paralogs
#echo $path
#echo $threads
#echo $grapple_dir
#echo "high and low set later in script"
#exit 0


# 1. Recreate PIRATE.all_alleles.tsv

if test -f ./PIRATE.all_alleles.tsv; then
	echo "PIRATE.all_alleles.tsv already exists"
else 
	mkdir ${grapple_dir}/
	$path/scripts/link_clusters_runner.pl -l ./loci_list.tab -t $thr_list -o ./${grapple_dir}/ -c ./co-ords/ -parallel $threads --all-alleles
	mv ${grapple_dir}/PIRATE.all_alleles.tsv ./PIRATE.all_alleles.tsv
fi


# 2. Create paralog files

if [[ $paralogs -eq 1 ]]; then
	# Recreate split_paralogs_loci.tab if necessary
	if test -f ./split_paralog_loci.tab; then
		echo "split_paralogs_loci.tab exists"
	else 
		$path/scripts/split_paralogs_runner.pl -p ./loci_paralog_categories.tab -l ./loci_list.tab -o ./ -t $threads
	fi
	
	wp="wp."
	
	if test -f ./PIRATE.all_alleles.wp.tsv; then
		echo "PIRATE.all_alleles.wp.tsv already exists"
	else
		mkdir ${grapple_dir}/
		${path}/scripts/link_clusters_runner.pl -l ./loci_list.tab -l ./split_paralog_loci.tab -t $thr_list -o ./${grapple_dir}/ -c ./co-ords/ --paralogs ./loci_paralog_categories.tab -e ./paralog_clusters.tab --parallel $threads --all-alleles
		mv ${grapple_dir}/PIRATE.all_alleles.tsv ./PIRATE.all_alleles.wp.tsv
        fi
else
	wp=""
fi

# 3. Split all_alleles by thresholds

if test -z $high; then
	# if not set, uses highest value to remove core genes (reduces data size for subsequent analyses) 
	high=$(cut -f7 ./PIRATE.all_alleles.${wp}tsv | sort -n | tail -1)
fi

if test -z $low; then
	# if not set, defaults to 1, to remove singleton genes (reduces data size for subsequent analyses)
	low=1
fi

for i in ${thr_list//,/ }
do
	awk -F"\t" -v thr=$i -v max=$high -v min=$low 'NR==1 {print $0}; $5==thr && $7<max && $7>min {print $0}' PIRATE.all_alleles.${wp}tsv > PIRATE.acc_alleles.${i}.${wp}genes_${low}-${high}.tsv 
done

# 4. Convert PIRATE files to binary format (necessary for pw_similarity.py)

for i in ${thr_list//,/ }
do
	$path/tools/convert_format/PIRATE_to_Rtab.pl -i ./PIRATE.acc_alleles.${i}.${wp}genes_${low}-${high}.tsv -o ./PIRATE.acc_alleles.${i}.${wp}genes_${low}-${high}.binary.tsv --low 0 --high 1 
done

