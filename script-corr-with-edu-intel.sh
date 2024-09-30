module load openblas ;module load fsl; export FSLOUTPUTTYPE=NIFTI2_GZ
fsl_glm -i /data/clusterfs/lag/users/sousoh/ukbb/genetic/img/t1/matrix-post-QC-abs-jac-iso-2/all-jac-absolute \
	-m  /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/t1-syn-template/qc-masks/stage_11_mask_roi_iso2 \
	-d /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/tbm/prs-matrices/dyslexia-lasso_col_283-edu-intel.mat \
	-c /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/dmri/contrasts/effect-pos-28-intel.con \
	--dat_norm  --out_cope=pheno-corr/intel_absolute-cope --out_varcb=pheno-corr/intel_absolute-varcope


fsl_glm -i /data/clusterfs/lag/users/sousoh/ukbb/genetic/img/t1/matrix-post-QC-abs-jac-iso-2/all-jac-absolute \
	-m  /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/t1-syn-template/qc-masks/stage_11_mask_roi_iso2 \
	-d /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/tbm/prs-matrices/dyslexia-lasso_col_283-edu-intel-headsize.mat \
	-c /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/dmri/contrasts/effect-pos-29-intel.con \
	--dat_norm  --out_cope=pheno-corr/intel_headsizenorm-cope --out_varcb=pheno-corr/intel_headsizenorm-varcope

fsl_glm -i /data/clusterfs/lag/users/sousoh/ukbb/genetic/img/t1/matrix-post-QC-abs-jac-iso-2/all-jac-absolute \
	-m  /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/t1-syn-template/qc-masks/stage_11_mask_roi_iso2 \
	-d /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/tbm/prs-matrices/dyslexia-lasso_col_283-edu-intel.mat \
	-c /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/dmri/contrasts/effect-pos-28-edu.con \
	--dat_norm  --out_cope=pheno-corr/edu_absolute-cope --out_varcb=pheno-corr/edu_absolute-varcope


fsl_glm -i /data/clusterfs/lag/users/sousoh/ukbb/genetic/img/t1/matrix-post-QC-abs-jac-iso-2/all-jac-absolute \
	-m  /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/t1-syn-template/qc-masks/stage_11_mask_roi_iso2 \
	-d /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/tbm/prs-matrices/dyslexia-lasso_col_283-edu-intel-headsize.mat \
	-c /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/dmri/contrasts/effect-pos-29-edu.con \
	--dat_norm  --out_cope=pheno-corr/edu_headsizenorm-cope --out_varcb=pheno-corr/edu_headsizenorm-varcope


