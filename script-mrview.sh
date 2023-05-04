mrview -voxelinfo 0 -noannotations -colourbar 1 \
	 -load /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/t1-syn-template/stage_11_syn_template_roi-t1-and-adc-composite.nii.gz -intensity_rang 0,3 \
	-overlay.load /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/tbm-prs-maps/dyslexia/randomise-pos_tstat1.nii.gz -overlay.colourmap 2 -overlay.threshold_min 3 -overlay.no_threshold_max -overlay.intensity 2,6 \
	 -overlay.load /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/tbm-prs-maps/dyslexia/randomise-neg_tstat1.nii.gz -overlay.colourmap 1 -overlay.threshold_min 3 -overlay.no_threshold_max -overlay.intensity 2,6 \	
	-overlay.load /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/tbm-prs-maps/dyslexia/randomise-pos_vox_corrp_tstat1_highres_minuslog10p.nii.gz -overlay.colour 0,1,0 -overlay.threshold_min 1.3 -overlay.no_threshold_max -overlay.intensity 0,4 \
	-overlay.load /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/tbm-prs-maps/dyslexia/randomise-neg_vox_corrp_tstat1_highres_minuslog10p.nii.gz -overlay.colour 186,85,255 -overlay.threshold_min 1.3 -overlay.no_threshold_max -overlay.intensity 0,4 \
	-fixel.load /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/fod-template/fixel-2-replica-with-maps/glm-z-lassosum-p195-background-fixel-visualization-mask.nii.gz \
 	-fixel.load /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/fod-template/fixel-2-replica-with-maps/randomise-dyslexia_positive_contrast_tstat1.nii.gz \
	-fixel.load /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/fod-template/fixel-2-replica-with-maps/randomise-dyslexia_negative_contrast_tstat1.nii.gz \
       -fixel.load /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/fod-template/fixel-2-replica-with-maps/randomise-pos_vox_corrp_tstat1_minuslog10p.nii.gz \
	-fixel.load /data/clusterfs/lag/users/sousoh/ukbb/pilot-1k/fod-template/fixel-2-replica-with-maps/randomise-neg_vox_corrp_tstat1_minuslog10p.nii.gz       
      
	       

