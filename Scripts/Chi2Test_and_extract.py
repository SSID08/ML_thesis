import pickle
import argparse
#import pandas as pd
import numpy as np
import os
import re
import subprocess as sp
from scipy.stats import chi2_contingency,false_discovery_control
from scipy.stats.contingency import odds_ratio

parser=argparse.ArgumentParser(description='Run chi square tests on shortlisted stop gained/lost variants')

parser.add_argument('--input',help="Pickle file with matrix",required=True)
parser.add_argument('--out_folder',help="output folder path",required=True)
#parser.add_argument('--out')
args=parser.parse_args()

file=pickle.load(open(args.input,"rb"))
DR_name=re.sub('.pkl','',os.path.basename(args.input))
matrix=np.array(file['Matrix'])
gene_locii=np.array(file['Column_IDs'])
pheno=np.array(file['Phenotype'])
p_vals=[]
odds_ratio_list=[]
#col_freqs=np.mean(matrix,axis=0)
#high_freq_variant_index=np.nonzero(col_freqs>0.1)
#high_freq_variant_freqs=col_freqs[high_freq_variant_index[0]]
#high_freq_variant_genes=gene_locii[high_freq_variant_index[0]]
positive_phenos=np.where(pheno==1)[0]
negative_phenos=np.where(pheno==0)[0]
positive_genos_matrix=matrix[positive_phenos,:]
negative_genos_matrix=matrix[negative_phenos,:]
positive_genos_sums=np.sum(positive_genos_matrix,axis=0)
negative_genos_sums=np.sum(negative_genos_matrix,axis=0)
#print(positive_genos_sums)
#print(negative_genos_sums)
#positive_genos_freq=np.mean(positive_genos_matrix,axis=0)
#negative_genos_freq=np.mean(negative_genos_matrix,axis=0)
#print(len(positive_genos_sums))
#print(len(positive_phenos))
#print(positive_genos_sums)
#print(positive_genos_sums[0])
#print(len(positive_phenos)-positive_genos_sums[0])
#print(len(positive_phenos)-positive_genos_sums[4])
# chi2_result=chi2_contingency(np.array([[0,2964],[1,7208]]))
# odds_r=odds_ratio(np.array([[0,2964],[1,7208]]))
# print(chi2_result)
# print(odds_r)
for i in range(0,len(positive_genos_sums)):
    variant_in_resistant=positive_genos_sums[i]+1
    no_variant_in_resistant=len(positive_phenos)-variant_in_resistant-1
    variant_in_susceptible=negative_genos_sums[i]+1
    no_variant_in_susceptible=len(negative_phenos)-variant_in_susceptible-1
    try:
        input=np.array([[variant_in_resistant,no_variant_in_resistant],[variant_in_susceptible,no_variant_in_susceptible]])
        chi2_result=chi2_contingency(input)
        odds_r=odds_ratio(input)
        odds_ratio_list.append(odds_r.statistic)
        p_vals.append(chi2_result.pvalue)
    except Exception as e: 
        print(e)

#f_name=os.path.join(args.out_folder,f'{DR_name}.bed')
#f=open(f_name,'w')
p_adj=false_discovery_control(ps=np.array(p_vals),method='bh')
regions_dict={}
for i in range(0,len(p_adj)):
    if p_adj[i]<0.05:
        regions_dict[str(gene_locii[i])]=[round(odds_ratio_list[i],3),"{:.3g}".format(p_adj[i])]
        #f.write(f'Chromosome\t{gene_locii[i]}\t{gene_locii[i]}\t{round(odds_ratio_list[i],3)}\t{"{:.3g}".format(p_adj[i])}\t{positive_genos_sums[i]}\t{negative_genos_sums[i]}\n')
regions_list=[f'Chromosome:{i}' for i in list(regions_dict.keys())]
#f.close()

out_dict={}

for l in sp.Popen(rf"bcftools view -r {','.join(regions_list)} ./Data/snpEff_annotated_vcf_DRonly.vcf.gz | bcftools query -f '%POS\t%AF\t%LOF\t%ANN\n'",shell=True,stdout=sp.PIPE).stdout:
    pos,AF,lof,ann=l.decode().strip().split()
    if lof!='.':
        try:
            pos=str(pos)
            OR=regions_dict[pos][0]
            p=regions_dict[pos][1]
            ann=ann.split(',')[0].split('|')
            var_type,gene,subs,protein_subs=ann[1],ann[3],ann[9],ann[10]
            if pos in out_dict:
                if AF>out_dict[pos][0]:
                    out_dict[pos]=[AF,var_type,gene,subs,protein_subs,OR,p]
            else:
                out_dict[pos]=[AF,var_type,gene,subs,protein_subs,OR,p]
        except Exception as e:
            print(f'Error {e} at position {pos}')

with open(f"{args.out_folder}/{DR_name}_Associated_variants.txt",'w') as f: 
    f.write("Position\tOdds_Ratio\tP_value\tGene\tVariant_Type\tNucleotide_Subs\tProtein_Subs\tAllele Frequency\n")
    for i in out_dict:
        try:
            arr=out_dict[i]
            f.write(f"{i}\t{arr[5]}\t{arr[6]}\t{arr[2]}\t{arr[1]}\t{arr[3]}\t{arr[4]}\t{arr[0]}\n")
        except Exception as e:
            print(f'Error {e} at position {i}')








# high_freq_in_pos=np.nonzero(positive_genos_freqs>0.1)
# df2=pickle.load(open("./Updated_Data/Pickle_files_for_LR/ethambutol.pkl","rb"))
# df3=pickle.load(open("./Data/Pickle_files_DRonly_positionOnly/ciprofloxacin.pkl","rb"))
# matrix=mat_file['Matrix']
# geneIDs=mat_file['GeneIDs']
# df=pickle.load(open("./Updated_Data/Transformed_dataframes_All_DR_SNPeff_high/pyrazinamide.pkl","rb"))['Matrix']
# df.shape
# cols=mat_file['Columns']
# pheno=mat_file['Phenotype']
# cols[0:10]
# genos=mat_file['geno_list']
# genos=np.arra
# len(mat_file['Columns'])
# matrix=np.array(matrix)
# pheno_file=mat_file['pheno_list']
# len(pheno_file)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          