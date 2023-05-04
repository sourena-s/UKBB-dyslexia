#!/bin/bash

index=$1
module load R
module load anaconda
conda activate ldsc
unset PYTHONPATH
unset PYTHONHOME

node_cores=$(echo "$(lscpu|awk '$1=="CPU(s):"{print $2}') - 4"|bc)
node_mem=$(free -g|awk 'NR==2{print $2}')
node_free_mem=$(free -g|awk 'NR==2{print $4}')

cond=()

cond[1]="dyslexia"
cond[2]="AD"
cond[3]="total-surf"
cond[4]="putamen"
cond[5]="SCZ"
cond[6]="ADHD"
cond[7]="ASD"
cond[8]="grade-E1"
cond[9]="grade-E2"
cond[10]="grade-E3"
cond[11]="grade-E4"
cond[12]="migraine"
cond[13]="PD"
cond[14]="BIP"
cond[15]="UKBB-GCSE"
cond[16]="UKBB-NVQ"
cond[17]="UKBB-PAIN"
cond[18]="UKBB-REACTION-TIME"
cond[19]="UKBB-VNR"
cond[20]="UKBB-VNR-comp1"
cond[21]="UKBB-VNR-comp2"
cond[22]="UKBB-VNR-comp3"
cond[23]="UKBB-VNR-comp4"
cond[24]="UKBB-VNR-comp5"
cond[25]="UKBB-VNR-comp6"
cond[26]="UKBB-VNR-comp7"
cond[27]="UKBB-VNR-comp8"
cond[28]="UKBB-VNR-comp9"
cond[29]="UKBB-VNR-comp10"
cond[30]="UKBB-VNR-comp11"
cond[31]="UKBB-VNR-comp12"
cond[32]="UKBB-VNR-comp13"
cond[33]="UKBB-hypertension"
cond[34]="UKBB-risk-taking"

#python=/usr/shared/apps/python/python-2.7.15/bin/python


echo "Running ${cond[index]} on a node with $((node_cores+4)) cores (reserving $node_cores cores) and $node_mem GB memory, $node_free_mem GB of which is free."
mkdir /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}

if [ ! -f  /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/ldsc_python_h2.log ];then

Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-sumstat-preproc.R $index
/home/sousoh/software/ldsc/munge_sumstats.py --sumstats /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/${cond[index]}/auto-sumstats.ma --out /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.munged --frq freq --merge-alleles /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/w_hm3.snplist

/home/sousoh/software/ldsc/ldsc.py --maf 0.01 --h2 /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/${cond[index]}/auto-sumstats.munged.sumstats.gz \
        --ref-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
        --w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
        --out /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/ldsc_python_h2

fi
