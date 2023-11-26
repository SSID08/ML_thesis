#! /bin/bash

cd ~/ML_thesis/Scripts/sh_scripts;

for file in ../../LOF_matrix/*;
do
base_name=$(basename "$file" .pkl);
echo $base_name
python ../Chi2Test_and_extract.py \
--input $file \
--anno_VCF ../../../VCF/DR_only_LOF.vcf.gz \
--out_folder ../../LOF_analysis_Results
done
