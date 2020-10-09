#!\bin\bash

set -oe pipefail

thr_list=$1
wp="wp."

for i in ${thr_list//,/ }
do
	files+=(PIRATE.acc_alleles.${i}.${wp}tsv)
	files+=(PIRATE.acc_alleles.${i}.${wp}binary.tsv)
done

echo "${files[@]}"

# for value in "${arr[@]}"
# do
#     echo $value
# done