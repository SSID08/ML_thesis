import subprocess as sp 
from collections import defaultdict
import numpy as np
import argparse
import pickle
import os
import re
import pandas as pd

parser=argparse.ArgumentParser(description='Get High impact variants genotype')
parser.add_argument('--bed_file',required=True)
parser.add_argument('--pheno',required=True)
parser.add_argument('--out_folder',required=True)
args=parser.parse_args()

tracker=defaultdict(dict)
gene_ids=[]
gene_names=[]

try: 
    with open(args.bed_file,'r') as f:
        for l in f:
            l_split=l.split()
            gene_ids.append(l_split[3])
            gene_names.append(l_split[4])
    print(gene_ids)
    print(gene_names)
except Exception as e: 
    print(f'Error {e} for file {args.bed_file}')

try:
    for l in sp.Popen(rf"bedtools intersect -a ./Data/snpEff_filterbyImpact.vcf.gz -wa -header -b {args.bed_file}" + r"| bcftools query -f '[%ANN\t%SAMPLE\t%GT\n]'",shell=True,stdout=sp.PIPE).stdout:
        ann,sample,gt=l.decode().strip().split()
        ann_split=ann.split(',')
        try:
            for i in ann_split:
                i_split=i.split('|')
                gene_name=i_split[3]
                gene_id=i_split[4]
                if (gene_name in gene_names) or (gene_id in gene_ids):
                    effect_anno=i_split[2]
                    key_name=str(f"{gene_name}_{effect_anno}")
                    if gt=="./." or gt=="0/0":
                        gt = 0
                    else:
                        gt = 1
                    try: 
                        tracker[sample]
                        try:
                            tracker[sample][key_name]+=gt
                        except KeyError:
                            tracker[sample][key_name]= 0 + gt
                    except KeyError:
                        tracker[sample][key_name]= 0 + gt
                    break
        except Exception as e:
            continue


    print(f'Genotype info done for file: {args.bed_file}')

    pheno={}

    with open(args.pheno,'r') as f: #Read phenotype file line by line
        next(f)
        for l in f:
            row=l.strip().split()
            if row[0] in tracker:
                pheno[row[0]]=int(row[1])

    print(f'Phenotype info done for file: {args.bed_file}')

    my_keys=list(tracker[list(tracker.keys())[0]].keys())
    out_list=[]
    Y=[]
    for s in pheno:
        #my_dict=tracker[s]
        out_list.append([tracker[s][key] for key in my_keys])
        Y.append(pheno[s])

    #Y = [pheno[s] for s in pheno]
    X=np.array(out_list)
    f_name=re.sub('.bed','',os.path.basename(args.bed_file))
    pickle.dump({'Matrix':X,'Phenotype':np.array(Y),'Column_IDs':my_keys},open(os.path.join(args.out_folder,f'{f_name}.pkl'),'wb'))

except Exception as e: 
    print(f'Exception : {e} for file {args.bed_file}')
    
    

        
