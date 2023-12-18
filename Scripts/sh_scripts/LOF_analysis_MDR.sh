#! /bin/bash

cd ~/ML_thesis/Scripts/sh_scripts;

##for drug in moxifloxacin ofloxacin levofloxacin linezolid;
# for file in ../../Per_Drug_bedFile/*;
# do
# base_name=$(basename "$file" .bed);
# echo $base_name
awk '{print $1}' ../../Sample_Phenotypes/XDR-TB_resistance.txt > ../../../tmp/XDR-TB.txt

python ../Extract_Loss_of_Function_variants.py \
--VCF_anno_file ../../../VCF/DR_only_LOF.vcf.gz \
--Sample_IDs ../../../tmp/XDR-TB.txt \
--pheno ../../Sample_Phenotypes/XDR-TB_resistance.txt \
--out_folder ../../LOF_matrix
#done
echo 'Done'
rm ../../../tmp/XDR-TB.txt
