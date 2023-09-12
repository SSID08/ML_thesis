import subprocess as sp 
#from collections import defaultdict
#import numpy as np
import argparse
#import pickle
#import os
#import re
#import pandas as pd

'''This script is a small helper program to merge Feature importance coefficients\
with SNPEff annotated information on the particular extracted variants of interest'''

parser=argparse.ArgumentParser(description='Inspect high-feature importance variants')
parser.add_argument('--vcf',required=True)
parser.add_argument('--bed_file',required=True)
parser.add_argument('--suffix',required=True)
parser.add_argument('--out_folder',required=True)
args=parser.parse_args()

out_dict={}
bed_positions={}
coefs=[]
count=0

for l in sp.Popen(rf"bedtools intersect -a {args.vcf} -wa -header -b {args.bed_file}" + r"| bcftools query -f '%AF\t%POS\t%ANN\n'",shell=True,stdout=sp.PIPE).stdout:
    count+=1
    AF,pos,ann=l.decode().strip().split()
    try:
        AF=float(AF)
    except:
        AF=float(0)
    pos=str(pos)
    try:
        ann=ann.split(',')[0].split('|')
        var_type,gene,subs,protein_subs=ann[1],ann[3],ann[9],ann[10]
        if pos in out_dict:
            if AF>out_dict[pos][0]:
                out_dict[pos]=[AF,var_type,gene,subs,protein_subs]
        else:
            out_dict[pos]=[AF,var_type,gene,subs,protein_subs]
    except Exception as e:
        # if pos in err_positions:
        #     if AF>err_positions[pos]:
        #         err_positions[pos]=AF
        # else:
        #     err_positions[pos]=AF
        print(f'Error {e} at position {pos}')


with open(args.bed_file,'r') as f: 
    for l in f:
        pos,coef=l.split()[1],l.split()[3]
        pos=str(pos)
        coef=round(float(coef),4)
        bed_positions[pos]=coef
        #bed_positions.append(pos)
        #coefs.append(coef)

with open(f"{args.out_folder}/{args.suffix}_feature_annotations.txt",'w') as f: 
    f.write("Coefficient\tPosition\tGene\tVariant_Type\tNucleotide_Subs\tProtein_Subs\tAlleleFrequency\n")
    for i in bed_positions:
        try:
            arr=out_dict[i]
            f.write(f"{bed_positions[i]}\t{i}\t{arr[2]}\t{arr[1]}\t{arr[3]}\t{arr[4]}\t{arr[0]}\n")
        except Exception as e: 
            continue
            #print(f'Problem with {bed_positions[i]} as {e}')
            
print(count)