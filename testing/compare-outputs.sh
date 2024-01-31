# returns 0 if the two arguments are the same
# returns 1 if the two arguments are different

OUTPUT1=$1
OUTPUT2=$2

SORTED1="output1_sorted.tsv"
SORTED2="output2_sorted.tsv"

sort ${OUTPUT1} > ${SORTED1}
sort ${OUTPUT2} > ${SORTED2}

exit $(($(diff ${SORTED1} ${SORTED2} -y --suppress-common-lines | wc -l) != 0))