import subprocess as sp 
from collections import defaultdict
import numpy as np
import pandas as pd
import argparse
import pickle
import os
import re

'''
This script creates feature matrix from VCF files for downstream ML analysis
'''
parser=argparse.ArgumentParser(description='Append MIC phenotype to MIC dataframe')
parser.add_argument('-MIC',help='MIC phenotype dataframe',required=True)
parser.add_argument('-pheno',help='Output folder',required=True)
args=parser.parse_args()

drug=(os.path.basename(args.pheno)).split('_')[0]
print(drug)
pheno_series={}

with open(args.pheno,'r') as f:
    for l in f: 
        sample,pheno=l.strip().split()
        pheno_series[sample]=pheno

pheno_series=pd.Series(pheno_series,name=drug,dtype='category')
print(pheno_series.cat.categories)
pheno_series=pheno_series.cat.reorder_categories(['NA','Sus','Low','High','Resistant'],ordered=True)
print(pheno_series.cat.categories)

df=pickle.load(open(args.MIC,'rb'))
print(f'Preivous shape: {df.shape}')
print(df.head())
out_df=pd.concat([df,pheno_series],axis=1,join='inner')
print(f'New shape: {out_df.shape}')
print(out_df.head())
#print(f'Sorted {out_df[drug].sort_values()}')

out_df.to_pickle(path=args.MIC)