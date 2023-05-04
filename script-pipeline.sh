#!/bin/bash

index=$1
module load R
module load anaconda
. ~/.bashrc

parallel="/home/sousoh/software/parallel/bin/parallel"
source /home/sousoh/software/parallel/bin/env_parallel.bash

node_cores=$(echo "$(lscpu|awk '$1=="CPU(s):"{print $2}') - 10"|bc)
node_mem=$(free -g|awk 'NR==2{print $2}')
node_free_mem=$(free -g|awk 'NR==2{print $4}')

conda deactivate

echo -n "Activating LDSC Anaconda environment.."
conda activate ldsc
echo "done."

unset PYTHONPATH
unset PYTHONHOME

cond=()

cond[1]="dyslexia"
cond[2]="AD-jansen-2019"
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
cond[35]="AD-kunkle"
cond[36]="handedness"
cond[37]="handedness-L"

#python=/usr/shared/apps/python/python-2.7.15/bin/python


echo "Running ${cond[index]} on a node with $((node_cores+4)) cores (reserving $node_cores cores) and $node_mem GB memory, $node_free_mem GB of which is free."
if [ ! -d /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]} ];then mkdir /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]} ;fi

Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-sumstat-preproc.R $index
/home/sousoh/software/ldsc/munge_sumstats.py --sumstats /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.ma --out /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.munged --frq freq --merge-alleles /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/w_hm3.snplist

/home/sousoh/software/ldsc/ldsc.py --maf 0.01 --h2 /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.munged.sumstats.gz \
        --ref-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
        --w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
        --out /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/ldsc_python_h2

bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-sbayesR.sh ${cond[index]} $node_cores
#did it fail?
if [ ! -f /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sbayesr.pgs.snpRes ];then
	cp /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/sbayesr.plainzero.snpRes /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sbayesr.pgs.snpRes
fi	

Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-ldpred2-calculate-hapmap.R $index $node_cores

###PRS-CS

conda deactivate
echo -n "activating PRS-CS Anaconda environment.."
conda activate prscs
echo "done."

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


if [ ! -d /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/prs-cs ];then mkdir /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/prs-cs; fi
cat /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.ma|cut -d" " -f1-3,5,7| \
awk '{if (NR==1) {print "SNP A1 A2 BETA P";}else{print $0;} }'| \
sed 's/ /\t/g' > /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.prs-cs

n_sample=$(cat /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.ma|awk 'NR>1{print $8}'| awk 'BEGIN{OFMT="%f" }{ sum += $1; n++ } END { if (n > 0) print int(sum / n); }')

phi=()
phi[1]="--seed=1"
phi[2]="--phi=1e-6"
phi[3]="--phi=1e-4"
phi[4]="--phi=1e-2"
phi[5]="--phi=1"

for i in {1..5};do
for chr in {1..22};do
echo -e "${phi[i]}\t$chr"
done
done | $parallel -j $node_cores --colsep "\t" python /home/sousoh/software/PRScs/PRScs.py --ref_dir /home/sousoh/software/PRScs/ldblk_ukbb_eur \
--bim_prefix=/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40 \
--sst_file=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/auto-sumstats.prs-cs \
--chrom={2} --out_dir=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/prs-cs/pgs {1} --n_gwas=$n_sample
#--n_burnin=1 --n_iter=2 \

for phi in auto 1e-06 1e-04 1e-02 1e+00;do
for chr in {1..22};do
cat /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/prs-cs/pgs_pst_eff_a1_b0.5_phi${phi}_chr${chr}.txt
done > /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/prs-cs/pgs-autosome-phi-${phi}.txt
done

rm /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/prs-cs/pgs_pst*txt
#### END PRS-CS


#Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-pgs-corr-all-ics.R ${cond[index]}
Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-pgs-corr-three-models.R ${cond[index]}

param_col_tbm=$(cat /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/lasso/optimisation_points.txt |awk 'NR==2{print $1}')
param_col_dmri=$(cat /data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/lasso/optimisation_points.txt |awk 'NR==2{print $2}')

echo "Creating TBM design matrix with SBayesR PGS"
echo "Creating dMRI design matrix with SBayesR PGS"

bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-prs-design-matrix-generator-multimodal-normalise.sh ${cond[index]} sbayesr tbm
bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-prs-design-matrix-generator-multimodal-normalise.sh ${cond[index]} sbayesr dmri

echo "Creating TBM design matrix with regularisation column $param_col_tbm"
echo "Creating dMRI design matrix with regularisation column $param_col_dmri"

bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-prs-design-matrix-generator-multimodal-normalise.sh ${cond[index]} $param_col_tbm tbm
bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-prs-design-matrix-generator-multimodal-normalise.sh ${cond[index]} $param_col_dmri dmri

randomise_script=/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/${cond[index]}/run_randomise.sh
echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
	-S /bin/bash -N rand_tbm_${cond[index]}_pos \
	/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} lasso_col_$param_col_tbm tbm pos" > $randomise_script

echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
	-S /bin/bash -N rand_tbm_${cond[index]}_neg \
	/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} lasso_col_$param_col_tbm tbm neg" >> $randomise_script


echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
        -S /bin/bash -N rand_dmri_${cond[index]}_pos \
        /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} lasso_col_$param_col_dmri dmri pos" >> $randomise_script

echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
        -S /bin/bash -N rand_dmri_${cond[index]}_neg \
        /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} lasso_col_$param_col_dmri dmri neg" >> $randomise_script

###SbayesR
echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
	-S /bin/bash -N rand_tbm_${cond[index]}_sbayesr_pos \
	/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} sbayesr tbm pos" >> $randomise_script

echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
	-S /bin/bash -N rand_tbm_${cond[index]}_sbayesr_neg \
	/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} sbayesr tbm neg" >> $randomise_script


echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
        -S /bin/bash -N rand_dmri_${cond[index]}_sbayesr_pos \
        /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} sbayesr dmri pos" >> $randomise_script

echo "qsub -q big.q -wd /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts \
        -S /bin/bash -N rand_dmri_${cond[index]}_sbayesr_neg \
        /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-LM-PRS.sh ${cond[index]} sbayesr dmri neg" >> $randomise_script



Rscript /data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/script-ggsave.R ${cond[index]}


