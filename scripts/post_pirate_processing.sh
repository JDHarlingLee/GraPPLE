#!/usr/bin/bash

# Author: JDHL

# Basic runner script for adapting the general PIRATE output for use in GraPPLE/Graphia
# Recreates some files (all_alleles.tsv), and then splits out files by identity threshold

# Can generate files with paralogs split or not (will generate relevant files if not included in the original PIRATE analysis)
# Files including split paralogs will be indicated with ".wp." in the file name

# Heavily relies on PIRATE's excellent adapter scripts, provided in the PIRATE repository
# Please see SionBayliss/PIRATE for more information on these

# read variables
while [ "$1" != "" ]; do
        case $1 in
                -t | --thr_list )	shift
					thr_list=$1
                                        ;;
	        --paralogs )		shift
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
		-m | --min_max )	shift
					min_max=$1
					;;
		-h | --help )           printf "\n------------------------------------------------\n"
					printf "\nThis is a script to adapt the standard output from PIRATE for use in GraPPLE/Graphia"
					printf "\nYou may prefer to run some of these actions individually, but this is a wrapper script for convenience\n"
					printf "\nEnsure you are executing this script within the output directory of your PIRATE pangenome analysis\n"
					printf "\n-t | threshold list (those used in PIRATE run, or a subset thereof)"
					printf "\n-d | output directory name - a directory where new files can be created to avoid overwriting originals. Default: GraPPLE"
					printf "\n-q | path to PIRATE scripts/ directory - necessary for using PIRATE adapter scripts. Default: Infer from first PIRATE executable in \$PATH environment variable."
					printf "\n-n | number of threads to use. Default: 2"
          				printf "\n-m | filter by the minimum and maximum number of genomes containing a gene (in the format minimum:maximum). Single values are interpreted as minimum. Default: 2"
					printf "\n--paralogs | flag to generate files with paralogs split. Default: off"
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

if test -z "$path"; then
	path=$(dirname $(dirname $(which PIRATE)))
	echo "PIRATE scripts directory location inferred: $path"
fi

if ! test -e "$path" || ! test -f $path/scripts/link_clusters_runner.pl || ! test -f $path/scripts/split_paralogs_runner.pl; then
        echo "ERROR: Path to PIRATE is invalid or is missing necessary adapter scripts"
        exit 1
fi

if test -f $path/bin/PIRATE_to_Rtab.pl; then
    Rtab_script_path="bin"
elif test -f $path/tools/convert_format/PIRATE_to_Rtab.pl; then
    Rtab_script_path="tools/convert_format"
else
    Rtab_script_path="None"
    echo "ERROR: Path to PIRATE is invalid or is missing necessary PIRATE_to_Rtab.pl adapter script"
    exit 1
fi

if test -z "$threads"; then
        threads=2
fi

if test -z "$grapple_dir"; then
	grapple_dir="GraPPLE"
fi

if test -d "$grapple_dir"; then
	echo "ERROR: Directory ${grapple_dir} already exists. Please choose a different directory name. Exiting to avoid overwriting previous data"
	exit 1
fi

 # Correct likely delimiter typos in the min_max parameter
min_max=$(echo ${min_max} | sed 's/[,;/-]/:/g')
if test -z "$min_max"; then
	min_max=2
elif test $(echo ${min_max} | fold -w1 | grep -vc "[0-9:]") -gt 0; then
        echo "ERROR: ${min_max} - min_max value may only contain numeric and delimiter characters. Accepted delimiters are : ; , - / characters."
	exit 1
fi

#For debugging:
#echo $thr_list
#echo $paralogs
#echo $path
#echo $threads
#echo $grapple_dir
#echo $min_max
#ls -lh $path/$Rtab_script_path/PIRATE_to_Rtab.pl
#exit 0

# 1. Create paralog files and PIRATE.all_alleles.tsv (if not already present)

if [[ $paralogs -eq 1 ]]; then
	# Create split_paralog_loci.tab if necessary
	if test -f ./split_paralog_loci.tab; then
		echo "split_paralog_loci.tab exists"
	else
		$path/scripts/split_paralogs_runner.pl -p ./loci_paralog_categories.tab -l ./loci_list.tab -o ./ -t $threads
	fi

	wp="wp."

	if test -f ./${grapple_dir}/PIRATE.all_alleles.wp.tsv; then
		echo "PIRATE.all_alleles.wp.tsv already exists"
	else
		mkdir ${grapple_dir}/
		${path}/scripts/link_clusters_runner.pl -l ./loci_list.tab -l ./split_paralog_loci.tab -t $thr_list -o ./${grapple_dir}/ -c ./co-ords/ --paralogs ./loci_paralog_categories.tab -e ./paralog_clusters.tab -parallel $threads --all-alleles
		mv ${grapple_dir}/PIRATE.all_alleles.tsv ${grapple_dir}/PIRATE.all_alleles.wp.tsv
        fi
else
	wp=""
	if test -f ./${grapple_dir}/PIRATE.all_alleles.tsv; then
		echo "PIRATE.all_alleles.tsv already exists"
	else
		mkdir ${grapple_dir}/
		$path/scripts/link_clusters_runner.pl -l ./loci_list.tab -t $thr_list -o ./${grapple_dir}/ -c ./co-ords/ -parallel $threads --all-alleles
		#mv ${grapple_dir}/PIRATE.all_alleles.tsv ./PIRATE.all_alleles.tsv
	fi
fi

# 2. Filter genes to reduce data size for subsequent analyses

 # calculates the number of delimiters in the --min_max value string
n_lowhigh_delimiters=$(echo ${min_max} | fold -w1 | grep -c ':')

# warns if the --min_max value contains more than two delimited values
if test ${n_lowhigh_delimiters} -gt 1 ; then
  echo "WARNING: More than two numerical values detected within --min_max parameter, using first two"
fi

# if set, extracts the minimum and maximum number of genes to include from before and after the delimiter
if test ${n_lowhigh_delimiters} -gt 0 ; then
  minimum=$(echo ${min_max} | cut -d: -f1)
  maximum=$(echo ${min_max} | cut -d: -f2)
  low=$(($minimum-1)) # calculates the highest integer below the minimum value
  high=$(($maximum+1)) # calculates the lowest integer above the min_maximum value
# if the min_max string contains illegal (non-numeric) characters, use default instead.
else
  minimum=$(echo ${min_max})
  # if not set, extracts maximum value from the PIRATE table
  maximum=$(cut -f7 ./${grapple_dir}/PIRATE.all_alleles.${wp}tsv | sort -n | tail -1)
  low=$(($minimum-1)) # calculates the highest integer below the minimum value
  high=$(($maximum+1)) # calculates the lowest integer above the min_maximum value
fi
echo " - filtering: gene clusters present in ${minimum}:${maximum} genomes will be retained"

# 3. Split all_alleles by thresholds

for i in ${thr_list//,/ }
do
	awk -F"\t" -v thr=$i -v max=$high -v min=$low 'NR==1 {print $0}; $5==thr && $7<max && $7>min {print $0}' ./${grapple_dir}/PIRATE.all_alleles.${wp}tsv > ./${grapple_dir}/PIRATE.acc_alleles.${i}.${wp}genes_${minimum}-${maximum}.tsv
done

# 4. Convert PIRATE files to binary format (necessary for pw_similarity.py)

for i in ${thr_list//,/ }
do
	$path/$Rtab_script_path/PIRATE_to_Rtab.pl -i ./${grapple_dir}/PIRATE.acc_alleles.${i}.${wp}genes_${minimum}-${maximum}.tsv -o ./${grapple_dir}/PIRATE.acc_alleles.${i}.${wp}genes_${minimum}-${maximum}.binary.tsv --low 0 --high 1
done
