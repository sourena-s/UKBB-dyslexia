gctb=/home/sousoh/software/gctb_2.03beta_Linux/gctb

cond=$1
n_cores=$2
input_ma=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/$cond/auto-sumstats.ma

if [ ! -d /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/$cond ]
then 
mkdir /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/$cond
fi

output_prs=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/$cond/auto-sbayesr.pgs

$gctb --sbayes R \
      --mldm /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/mldm-list.txt \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $input_ma \
     --chain-length 50000 \
     --burn-in 10000 \
     --out-freq 10 \
     --out $output_prs \
     --p-value 1 \
     --unscale-genotype \
     --exclude-mhc \
     --thread $n_cores \
     --no-mcmc-bin \
     --rsq 0.9


