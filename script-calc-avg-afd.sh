#!/bin/bash
module load openblas
module load fsl
export FSLOUTPUTTYPE=NIFTI2_GZ

m_path=/data/clusterfs/lag/users/sousoh/ukbb/dmri-all/all-subjects-symlinks
subject=$1

fslstats ${m_path}/${subject}/afd.nii.gz -m > ${m_path}/${subject}/avg_afd.txt
