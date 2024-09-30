#$ -S /bin/bash
module load anaconda
module load openblas
source activate /home/sousoh/conda-env/mrtrix



for i in eddy-nonrotated-subs eddy-rotated-subs rescan;do

      cat ${i}-qced.txt|while read s;do

if [ -f $i/$s/mrtrix/fixel_in_template/afd.nii.gz ];then echo "skipping $i/$s"
else

      echo "Submitting $i/$s"

      echo " echo \"Analysis started at: \$(date) \" > ${i}/${s}/mrtrix/log.txt" > ${i}/${s}/mrtrix/nifti-script.sh
  echo " module load anaconda" >> ${i}/${s}/mrtrix/nifti-script.sh
  echo " module load openblas" >> ${i}/${s}/mrtrix/nifti-script.sh
  echo " source activate /home/sousoh/conda-env/mrtrix" >> ${i}/${s}/mrtrix/nifti-script.sh
  echo " mrconvert -info $i/$s/mrtrix/fixel_in_template/afd.mif $i/$s/mrtrix/fixel_in_template/afd.nii.gz" >> ${i}/${s}/mrtrix/nifti-script.sh

#  bash ${i}/${s}/mrtrix/nifti-script.sh
qsub -q single15.q -wd $(pwd) -o ${i}/${s}/mrtrix -e ${i}/${s}/mrtrix -S /bin/bash ${i}/${s}/mrtrix/nifti-script.sh



fi


done

done


