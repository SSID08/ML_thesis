import numpy as np
import pandas as pd
import pickle
import os
import argparse

parser=argparse.ArgumentParser(description='Write MIC DF to CSV')
parser.add_argument('-MIC',help='MIC phenotype dataframe',required=True)
parser.add_argument('-drug',help='Output folder',required=True)
args=parser.parse_args()

df=pickle.load(open(args.MIC,'rb'))

