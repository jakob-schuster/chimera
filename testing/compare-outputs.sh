# returns 0 if the two arguments are the same
# returns 1 if the two arguments are different

OUTPUT1=$1
OUTPUT2=$2



exit $(($(diff ${OUTPUT1} ${OUTPUT2} -y --suppress-common-lines | wc -l) != 0))