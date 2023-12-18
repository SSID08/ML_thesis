import os 
import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import re
import argparse

parser=argparse.ArgumentParser(description='Plot feature importances of model')

parser.add_argument('--model',help="Input model",required=True)
parser.add_argument('--df',help='Dataframe of variant counts',required=True)
parser.add_argument('--out_folder',help="output folder path",required=True)

args=parser.parse_args()
DR_name=re.sub('.pkl','',os.path.basename(args.model))

importances=pickle.load(open(args.model,'rb'))['Model'].feature_importances_
positions=pickle.load(open(args.df,'rb'))['Columns']
sorted_indices=importances.argsort()[::-1][0:10]
importances=importances[sorted_indices]
positions=positions[sorted_indices]

for l in sp.open(r'')



