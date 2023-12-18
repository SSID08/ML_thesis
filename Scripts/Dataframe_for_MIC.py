import subprocess as sp 
from collections import defaultdict
import numpy as np
import pandas as pd
import argparse
import os
#import re

'''
Create Dataframes for MIC analysis 
'''
parser=argparse.ArgumentParser(description='Create feature dataframe for MIC analysis')
parser.add_argument('--bed',required=True)
parser.add_argument('--VCF',help='Input VCF file',required=True)
parser.add_argument('--pheno',help='Output folder',required=True)

args=parser.parse_args()

genos=defaultdict(dict)
drug=os.path.basename(args.bed).split('.')[0]

for l in sp.Popen(rf"bedtools intersect -a {args.VCF} -b {args.bed} -wa -header | \
                  bcftools query -i '(TYPE=" + r'"snp" | TYPE="indel")' + r"'" + f" -f '[%POS\t%SAMPLE\t%GT\n]'",shell=True,stdout=sp.PIPE).stdout:
    pos,sample,gt=l.decode().strip().split() # split each line on tab
    pos=int(pos)
    if gt=="./." or gt=="0/0":#infer genotype value from line
        gt = 0
    else:
        gt = 1
    try:
        #See if sample key exists in dictionary 
        genos[sample]
        try:
            #See if 'pos' key exists for that particular sample
            if genos[sample][pos]==0 and gt==1: # only if curent genotype from line is 1 and existing call in dictionary is zero
                # add the genotype to the dictionary
                genos[sample][pos] = gt
        except KeyError: # if 'pos' key not found in dictionary 
            # add the genotype to the dictionary
            genos[sample][pos] = gt 
    except KeyError: # if 'sample' key not found in dictionary
        genos[sample][pos] = gt             
        # add the genotype to the dictionary    

print('Genotype data parsed')

try: 
    samples=list(genos.keys())
    genome_positions=list(genos[list(genos.keys())[0]].keys())
    X=[]
    for s in samples:
        X.append([genos[s][key] for key in genome_positions])

except Exception as e:
    print(e)

try:
    genos_matrix=np.array(X)
    print(genos_matrix[0:5,0:5])
    pd_df=pd.DataFrame(data=genos_matrix,columns=genome_positions,index=np.array(samples))
    #pd_df.to_pickle(path=os.path.join('~/ML_thesis','MIC_dfs',''))
except Exception as e: 
    print(f'Exception {e} occured while building dataframe')

pheno_series={}
try:
    with open(args.pheno,'r') as f:
        for l in f: 
            sample,pheno=l.strip().split()
            pheno_series[sample]=pheno

    pheno_series=pd.Series(pheno_series,name=drug,dtype='category')
    print(pheno_series.cat.categories)
    print(f'Previous shape {pd_df.shape}')
    pd_df=pd.concat([pd_df,pheno_series],axis=1,join='inner')
    print(f'New shape {pd_df.shape}')

except Exception as e:
    print(f'Exception {e} occured during phenotype addition for drug {drug}')

pd_df.to_csv(os.path.join('~/ML_thesis','MIC_dfs',f'{drug}_MIC_df.csv'))
