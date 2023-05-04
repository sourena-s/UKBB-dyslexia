#!/bin/bash
i=$1
cond_name=$2

module load R
module load anaconda
conda activate ldsc
unset PYTHONPATH
unset PYTHONHOME

#Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-sumstat-preproc.R $i

#python=/usr/shared/apps/python/python-2.7.15/bin/python

/home/sousoh/software/ldsc/munge_sumstats.py --sumstats /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/$cond_name/auto-sumstats.ma \
        --out /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/$cond_name/auto-sumstats.munged --frq freq --merge-alleles /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/w_hm3.snplist

/home/sousoh/software/ldsc/ldsc.py --maf 0.01 --h2 /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/$cond_name/auto-sumstats.munged.sumstats.gz \
        --ref-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
        --w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
        --out /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/$cond_name/ldsc_python_h2

