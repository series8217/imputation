#
# tests randomly perturbed sample haplotypes then runs an imputation
#

infile=$1
hapfile=$2
legfile=$3

cat $infile | tr " " "\n" | ./randperturb.sh per_$infile 5

./imputation per_$infile $hapfile $legfile imp_$infile readable_imp_$infile".txt"

cmp -l $infile imp_$infile | wc -l
echo " errors"
