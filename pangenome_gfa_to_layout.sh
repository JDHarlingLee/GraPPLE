#!\bin\bash

DATETIME="`date '+%Y%m%d%H%M%S'`"

while [ "$1" != "" ]; do
	case $1 in
		-pg | --pangenome_gfa )	shift
					pangenome_gfa=$1
					;;
		-m | --metadata ) 	shift
					metadata=$1
					;;
		-e | --edges)		shift
					edges=$1
					;;
		-t | --tidy )		tidy=1
					;;
		-h | --help )		echo "Ensure that all files are in the working directory"
					echo "-pg pangenome_gfa | -e edges | -m metadata | -t tidy mode | -h help"
					exit
					;;
		* )			echo "Use -h | --help"
					exit 1
	esac
	shift
done

grep "^S" $pangenome_gfa > nodeclass_${DATETIME}.tmp

# nodeclass.layout
# sed removes RC:i: from in front of number of samples (####TODO clarify what this means!)
sed 's/RC:i://g' nodeclass_${DATETIME}.tmp | awk -v q="\"" {'print "//NODECLASS""\t"q$2q"\t"q$4q"\t"q"number_samples"q'} > nodeclass_${DATETIME}.layout

# edges.layout
# sed removes - from node; current graph is unable to display directionality of synteny
sed 's/-//g' $edges | awk -v q="\"" {'print q$1q"\t"q$2q"\t"q$3q'} > edges_${DATETIME}.layout

# metadata
awk -F"\t" -v q="\"" {'print "//NODECLASS""\t"q$1q"\t"q$3q"\t"q"consensus_gene_name"q'} $metadata > metadata_${DATETIME}.layout
awk -F"\t" -v q="\"" {'print "//NODECLASS""\t"q$1q"\t"q$4q"\t"q"consensus_product"q'} $metadata >> metadata_${DATETIME}.layout
awk -F"\t" -v q="\"" {'print "//NODECLASS""\t"q$1q"\t"q$7q"\t"q"number_genomes"q'} $metadata >> metadata_${DATETIME}.layout
awk -F"\t" -v q="\"" {'print "//NODECLASS""\t"q$1q"\t"q$16q"\t"q"all_products"q'} $metadata >> metadata_${DATETIME}.layout
awk -F"\t" -v q="\"" {'print "//NODECLASS""\t"q$1q"\t"q$17q"\t"q"all_gene_names"q'} $metadata >> metadata_${DATETIME}.layout

# print layout file
cat edges_${DATETIME}.layout nodeclass_${DATETIME}.layout metadata_${DATETIME}.layout > pangenome_${DATETIME}.layout

echo "-----------------------------"
echo ""
echo "Pangenome Layout file printed"
echo ""
echo "-----------------------------"
echo ""
echo ""

rm nodeclass_${DATETIME}.tmp

if [ "$tidy" = "1" ]; then
	rm nodeclass_${DATETIME}.layout
	rm edges_${DATETIME}.layout
	rm metadata_${DATETIME}.layout
	echo "------------------------------------------------"
	echo ""
	echo "Tidy mode active: removed all intermediate files"
	echo ""
	echo "------------------------------------------------"
fi

