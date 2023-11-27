import subprocess as sp 
from collections import defaultdict
import numpy as np
import argparse
import pickle
import os
import re
#import pandas as pd

'''
Extract loss of function variants
'''
parser=argparse.ArgumentParser(description='Extract loss of function variants')
parser.add_argument('--VCF_anno_file',required=True)
parser.add_argument('--Sample_IDs',required=True)
parser.add_argument('--bed_file',required=False)
parser.add_argument('--pheno',required=True)
parser.add_argument('--out_folder',required=True)
args=parser.parse_args()

tracker=defaultdict(dict)
var_names=[]

# try: 
#     with open(args.bed_file,'r') as f:
#         for l in f:
#             l_split=l.split()
#             gene_ids.append(l_split[3])
#             gene_names.append(l_split[4])
#     print(gene_ids)
#     print(gene_names)
# except Exception as e: 
#     print(f'Error {e} for file {args.bed_file}')

if args.bed_file:
    cmd= f"bedtools intersect -a {args.VCF_anno_file} -wa -header -b {args.bed_file}" + f"| bcftools query -S {args.Sample_IDs} --force-samples -f '[%POS\t%LOF\t%SAMPLE\t%GT\n]'"
else: 
    cmd= f"bcftools query -S {args.Sample_IDs} --force-samples -f '[%POS\t%LOF\t%SAMPLE\t%GT\n]' {args.VCF_anno_file}"

try:
    for l in sp.Popen(cmd,shell=True,stdout=sp.PIPE).stdout:
        pos,lof,sample,gt=l.decode().strip().split()
        pos=int(pos)
        if lof != '.':
            if gt=="./." or gt=="0/0":
                gt = 0
            else:
                gt = 1
            try: 
                tracker[sample]
                try:
                    tracker[sample][pos]=gt
                except KeyError:
                    tracker[sample][pos]= gt
            except KeyError:
                tracker[sample][pos]= gt
    
    print(len(tracker))

    print(f'Genotype info done for file: {args.pheno}')

    pheno={}

    with open(args.pheno,'r') as f: #Read phenotype file line by line
        next(f)
        for l in f:
            row=l.strip().split()
            if row[0] in tracker:
                pheno[row[0]]=int(row[1])

    print(f'Phenotype info done for file: {args.pheno}')

    my_keys=list(tracker[list(tracker.keys())[0]].keys())
    out_list=[]
    Y=[]
    for s in pheno:
        #my_dict=tracker[s]
        out_list.append([tracker[s][key] for key in my_keys])
        Y.append(pheno[s])

    #Y = [pheno[s] for s in pheno]
    X=np.array(out_list)
    f_name=re.sub('.txt','',os.path.basename(args.pheno))
    pickle.dump({'Matrix':X,'Phenotype':np.array(Y),'Column_IDs':my_keys},open(os.path.join(args.out_folder,f'{f_name}.pkl'),'wb'))

except Exception as e: 
    print(f'Exception : {e} for file {args.pheno}')
    
    

        
