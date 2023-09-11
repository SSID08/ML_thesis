import pickle
import argparse
import numpy as np
import pandas as pd
import os
import re
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split

parser=argparse.ArgumentParser(description='Run Permutation Feature Importance')

parser.add_argument('--input',help="pickle input matrix and classification labels",required=True)
parser.add_argument('--out_folder',help="output folder path",required=True)
parser.add_argument('--feature_set',help="identifier for your input feature set",required=True)
parser.add_argument('--method',required=True)

args=parser.parse_args()
DR_name=re.sub('.pkl','',os.path.basename(args.input))

items=pickle.load(open(args.input,"rb"))
X=np.array(items["Matrix"])
Y=np.array(items["Phenotype"])
feature_names=np.array(items['Columns'])
#weight=((Y==0).sum())/(Y.sum())

#print(X.shape)
#print(Y.shape)
#min_samples_split=0.005,min_samples_leaf=0.005
rf=RandomForestClassifier(n_estimators=100,max_features=0.5,max_samples=0.7,min_samples_split=0.0001,min_samples_leaf=0.0001,class_weight='balanced')
X_train,X_test,Y_train,Y_test=train_test_split(X,Y,test_size=0.2,stratify=Y,shuffle=True)
rf=rf.fit(X_train,Y_train)
r = permutation_importance(rf, X_test, Y_test,n_repeats=20,random_state=12,scoring='roc_auc')
print('Permutation done')
with open(os.path.join(args.out_folder,f'{args.method}_{DR_name}_{args.feature_set}_permutationImportance.bed'),'w') as f:
    for i in r.importances_mean.argsort()[::-1]:
        if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
            name=str(feature_names[i])
            if name.startswith('pass'):
                position=re.sub('.*__','',name)
                f.write(f'Chromosome\t{position}\t{position}\t{r.importances_mean[i]}\t{r.importances_std[i]}\n')
                #f.write(f'{feature_names[i]}\t{r.importances_mean[i]}\t{r.importances_std[i]}\n')
            else:
                lineage=re.sub('.*Lineage_','',name)
                f.write(f'Lineage{lineage}\t1\t1\t{r.importances_mean[i]}\t{r.importances_std[i]}\n')

                