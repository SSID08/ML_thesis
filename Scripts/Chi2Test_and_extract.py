import pickle
import argparse
#import pandas as pd
import numpy as np
import os
import re
import subprocess as sp
from scipy.stats import fisher_exact,false_discovery_control
from scipy.stats.contingency import odds_ratio
import sys

parser=argparse.ArgumentParser(description='Run chi square tests on shortlisted stop gained/lost variants')

parser.add_argument('--input',help="Pickle file with matrix",required=True)
parser.add_argument('--out_folder',help="output folder path",required=True)
parser.add_argument('--anno_VCF',help='SNPeff annotated VCF',required=True)
#parser.add_argument('--out')
args=parser.parse_args()

file=pickle.load(open(args.input,"rb"))
DR_name=re.sub('.pkl','',os.path.basename(args.input))
matrix=np.array(file['Matrix'])
gene_locii=np.array(file['Column_IDs'])
pheno=np.array(file['Phenotype'])
p_vals=[]
odds_ratio_list=[]
cont_matrices=[]
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

num_pos=len(positive_phenos)
num_neg=len(negative_phenos)
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
locii_inc=[]
for i in range(0,len(positive_genos_sums)):
    variant_in_resistant=positive_genos_sums[i]+1
    no_variant_in_resistant=num_pos-variant_in_resistant+2
    variant_in_susceptible=negative_genos_sums[i]+1
    no_variant_in_susceptible=num_neg-variant_in_susceptible+2
    if variant_in_resistant>4:
        try:
            input=np.array([[variant_in_resistant,no_variant_in_resistant],[variant_in_susceptible,no_variant_in_susceptible]])
            fisher_result=fisher_exact(input)
            odds_r=odds_ratio(input)
            odds_ratio_list.append(odds_r.statistic)
            p_vals.append(fisher_result.pvalue)
            cont_matrices.append(input)
            locii_inc.append(i)
        except Exception as e: 
            print(e)

#f_name=os.path.join(args.out_folder,f'{DR_name}.bed')
#f=open(f_name,'w')
p_adj=false_discovery_control(ps=np.array(p_vals),method='bh')
regions_dict={}
for i in range(0,len(p_adj)):
    if p_adj[i]<0.05:
        regions_dict[str(gene_locii[locii_inc[i]])]=[round(odds_ratio_list[i],3),"{:.3g}".format(p_adj[i]),cont_matrices[i]]
        #f.write(f'Chromosome\t{gene_locii[i]}\t{gene_locii[i]}\t{round(odds_ratio_list[i],3)}\t{"{:.3g}".format(p_adj[i])}\t{positive_genos_sums[i]}\t{negative_genos_sums[i]}\n')
#regions_list=[f'Chromosome:{i}' for i in list(regions_dict.keys())]
#f.close()

if not (len(odds_ratio_list)==len(locii_inc)==len(p_adj)):
    print('Something has gone wrong with the indices')
    sys.exit(0)



# print(f'Number of gene Locii: {len(gene_locii)}')
# print(f'Number if OR calculated : {len(odds_ratio_list)}')
# print(f'Number of locii included: {len(locii_inc)}')
# print(f'Number of p-vals calculated : {len(p_adj)} ')
# print(regions_dict)

out_dict={}
non_sig_vals=0
for l in sp.Popen(rf"bcftools query -f '%POS\t%AF\t%LOF\t%ANN\n' {args.anno_VCF}",shell=True,stdout=sp.PIPE).stdout:
    pos,AF,lof,ann=l.decode().strip().split()
    AF=float(AF)
    if lof!='.':
        try:
            pos=str(pos)
            OR=regions_dict[pos][0]
            p=regions_dict[pos][1]
            matrix=regions_dict[pos][2]
            ann=ann.split(',')[0].split('|')
            var_type,gene,subs,protein_subs=ann[1],ann[3],ann[9],ann[10]
            if pos in out_dict:
                if AF>out_dict[pos][0]:
                    out_dict[pos]=[AF,var_type,gene,subs,protein_subs,OR,p,lof,matrix]
            else:
                out_dict[pos]=[AF,var_type,gene,subs,protein_subs,OR,p,lof,matrix]
        except Exception as e:
            non_sig_vals+=1

with open(f"{args.out_folder}/{DR_name}_Associated_variants.txt",'w') as f: 
    f.write("Position\tOdds_Ratio\tMatrix\tP_value\tGene\tVariant_Type\tNucleotide_Subs\tProtein_Subs\tAllele Frequency\tLOF\n")
    for i in out_dict:
        try:
            arr=out_dict[i]
            f.write(f"{i}\t{arr[5]}\t{np.array2string(arr[8].flatten(),separator=',')}\t{arr[6]}\t{arr[2]}\t{arr[1]}\t{arr[3]}\t{arr[4]}\t{str(arr[0])}\t{arr[7]}\n")
        except Exception as e:
            print(f'Error {e} at position {i}')

print(f'Done for drug {DR_name}. Number of non sig values is {non_sig_vals}')




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
# len(pheno_file)