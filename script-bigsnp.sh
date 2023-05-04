#!/bin/bash
module load R
varlist="/data/clusterfs/lag/users/sousoh/ukbb/genetic/UKBB-clump/var-clump-1e-3.txt"

Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/script-bigsnp.R

awk '{print $1}' $varlist| while read var;do 
bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-variant-design-matrix-generator.sh $var
done
