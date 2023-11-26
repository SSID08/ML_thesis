#! /bin/bash

cd ~/ML_thesis/Scripts/sh_scripts;

##for drug in moxifloxacin ofloxacin levofloxacin linezolid;
for file in ../../Per_Drug_bedFile/*;
do
base_name=$(basename "$file" .bed);
echo $base_name
python ../Extract_Loss_of_Function_variants.py \
--VCF_anno_file ../../../VCF/DR_only_LOF.vcf.gz \
--Sample_IDs ../../Sample_info/${base_name}Samples.txt \
--bed_file ../../Per_Drug_bedFile/${base_name}.bed \
--pheno ../../Sample_Phenotypes/${base_name}Sample_Phenotypes.txt \
--out_folder ../../LOF_matrix
done
