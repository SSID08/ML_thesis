#! /bin/bash
cd ~/ML_thesis/Scripts/sh_scripts;

awk '{print $1}' ../../MIC_Phenotypes/INH_categorical_MIC_phenotypes.txt > ../../../tmp/MIC-pheno.txt

bcftools view -S ../../../tmp/MIC-pheno.txt --force-samples -Ou ../../../VCF/snpEff_annotated_vcf_DRonly.vcf.gz | \
bcftools view -q 0.005:nonmajor -Oz -o ../../../VCF/DRonly_filtered_0.005.vcf.gz

rm ../../../tmp/MIC-pheno.txt


