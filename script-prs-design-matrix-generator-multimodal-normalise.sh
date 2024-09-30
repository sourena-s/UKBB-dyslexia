#!/bin/bash


condition=$1 # e.g. dyslexia
regularisation_col=$2
modality=$3 #  "tbm" or "dmri"

design_base=/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/$modality

if [ "$regularisation_col" == "sbayesr" ];then
prsname="sbayesr"
output_f=/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/$modality/prs-matrices/$condition-$prsname
column_data=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${condition}/temp-sbayesr_pgs-col-$modality
paste -d" " /data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt \
	"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${condition}/auto-sbayesr.pgs" |cut -d" " -f1,3 |sort -k1,1 > $column_data
join -j1 ${design_base}/c0-post-QC $column_data |awk '{print $3}' > ${output_f}-col

elif [ "$regularisation_col" == "cs_auto" ];then

prsname="cs_auto"
output_f=/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/$modality/prs-matrices/$condition-$prsname
column_data=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${condition}/temp-cs_auto_pgs-col-$modality
paste -d" " /data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt \
	"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${condition}/prs-cs/subj_pgs_phi_auto.pgs" |cut -d" " -f1,3 |sort -k1,1 > $column_data
join -j1 ${design_base}/c0-post-QC $column_data |awk '{print $3}' > ${output_f}-col

else
prsname=lasso_col_$regularisation_col # e.g. 0.001

#PRS threshs: 0.001 0.01 0.1 10 15 1 20 25 30 35 40 45 50 55 5 60 65 70 75 80 85 90 95 100
#threshold can be pgs name, e.g. lassosum-par195
output_f=/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/$modality/prs-matrices/$condition-$prsname
#column_data=/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/prs-thresholds/${condition}/${threshold}.txt.sorted

#a<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/dyslexia/component-results/master_df.txt',header=T)
#b<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/dyslexia/lasso/grid-hapmap.params',sep=",",header=T)

cut -d, -f$regularisation_col /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${condition}/lasso/lasso-hapmap.pgs > temp-pgs

column_data=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${condition}/temp-lasso_pgs-col-$modality
paste -d" " /data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt temp-pgs |cut -d" " -f1,3 |sort -k1,1 > $column_data
#paste -d" " /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/AD/ldpred2/subject_id temp-pgs|sort -k1,1 > temp-pgs-col

join -j1 ${design_base}/c0-post-QC $column_data |awk '{print $3}' > ${output_f}-col
fi


if [ "$modality" == "tbm" ];then

paste -d" " ${design_base}/c0-post-QC ${output_f}-col \
	${design_base}/c0-5-coded ${design_base}/c7-pca_ten \
	${design_base}/design-agesq-ageXsex-agesqXsex ${design_base}/base-genoBatch |awk '$2=="OK"{print $0" 1"}' |cut -d" " -f3,5-  \
> ${output_f}.txt

#flirt
#sex
#age
#mri_site
#scan-date
#pca_10
#age_sq
#age_X_sex
#age_sq_X_sex
#geno_batch
#intercept

elif [ "$modality" == "dmri" ];then

paste -d" " ${design_base}/c0-post-QC ${output_f}-col \
	${design_base}/c1-sex-age-site-sitedate-rescanDateAgecorrected ${design_base}/c2-genomic-PCs \
   	${design_base}/c3-AgeSq-SexByAge-SexByAgeSq ${design_base}/c4-genoarray ${design_base}/c5-inHouseEddy \
	|sed 's/,/ /g' |awk '$2=="OK"'|cut -d" " -f3-|awk '{print $0" 1"}' > ${output_f}.txt


else
	echo "unknown modality! Options are tbm or dmri"

fi

Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/normalise-design-mat.R ${output_f}.txt
Text2Vest ${output_f}.txt.norm ${output_f}.mat

rm ${output_f}-col ${output_f}.txt


