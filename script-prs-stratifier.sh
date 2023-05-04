input_ma=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia-sumstats-clean-2.ma
input_prs=/data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/prs-dyslexia.snpRes
k=$1


if [ ! -f ${input_ma}.sorted ];then
awk 'NR>1' $input_ma | sort -k1,1 > ${input_ma}.sorted
fi

if [ ! -f  ${input_prs}.sorted ];then
awk 'NR>1' $input_prs | sort -k2,2 > ${input_prs}.sorted
fi


n_sign_concord=$(join -1 1 -2 2 ${input_ma}.sorted ${input_prs}.sorted|awk '$15 * $5 > 0'|wc -l)

n_tot_snp=$(cat $input_prs|wc -l)

echo "as.integer(${n_sign_concord} * ${k} / 100)" > ${input_prs}.${k}.R
n_top_snp=$(Rscript ${input_prs}.${k}.R |awk '{print $2}')

rm ${input_prs}.${k}.R

echo "number of total snp: $n_tot_snp"
echo "number of sign concording snp: $n_sign_concord"
echo "at $k percentile, $n_top_snp will be kept."

if [ ! -d ${input_prs}.stratified ];then
mkdir ${input_prs}.stratified
fi

echo "Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP LastSampleEff" > ${input_prs}.stratified/top-${k}-percentile
join -1 1 -2 2 ${input_ma}.sorted ${input_prs}.sorted| awk 'function abs(v) {return v < 0 ? -v : v} {print abs($15)" "$0}'| awk '$16 * $6 > 0'|sort -k1,1 -g -r |cut -d" " -f2,10- |head -n $n_top_snp >> ${input_prs}.stratified/top-${k}-percentile

#rm ${input_ma}.sorted ${input_prs}.sorted
