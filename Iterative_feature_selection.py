import pickle
#from collections import defaultdict
import argparse
import xgboost as xgb
#from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score,recall_score,make_scorer,roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_validate
#from sklearn.metrics import precision_score,recall_score
#from sklearn.metrics import confusion_matrix
import numpy as np
import os
import re
#import json

#genos=defaultdict(dict)

parser=argparse.ArgumentParser(description='Run ML training with iterative downstream feature selection')

parser.add_argument('--input',help="pickle input matrix and classification labels",required=True)
parser.add_argument('--out_folder',help="output folder path",required=True)
#parser.add_argument('--scale_pos_weight',help='Scale pos weight parameter',required=False)
#parser.add_argument('--feature_set',help="identifier for your input feature set",required=True)
parser.add_argument('--max_iter',help='Maximum number of iterations',required=True)

args=parser.parse_args()
DR_name=re.sub('.pkl','',os.path.basename(args.input))
max_iter=int(args.max_iter)
items=pickle.load(open(args.input,"rb"))
X=items["Matrix"]
Y=items["Phenotype"]
#weight=((Y==0).sum())/(Y.sum())
#base_score=Y.sum()/len(Y)
cols=items['Columns']
specificity_scorer=make_scorer(recall_score,pos_label=0)
NPV_scorer=make_scorer(precision_score,pos_label=0)
scorer={'recall':'recall','precision':'precision','specificity':specificity_scorer,'NPV':NPV_scorer,'accuracy':'accuracy','roc_auc':'roc_auc'}
xgb_classifier=xgb.XGBClassifier(objective='binary:logistic',booster='gbtree',learning_rate=.1)
X_train,X_test,Y_train,y_test=train_test_split(X,Y,test_size=0.3,stratify=Y,shuffle=True)
n_features=len(cols)
iter=0
best_roc,best_sens,best_spec,best_roc_n,best_sens_n,best_spec_n=[0]*6
best_roc_feature_set=[]
f=open(os.path.join(args.out_folder,DR_name,f"{DR_name}_Feature_Selection_Results"+".txt"),'w')
while n_features>100 and iter<max_iter:
    f.write(f'{X_train.shape}\n')
    base_model=xgb_classifier.fit(X_train,Y_train)
    pred=base_model.predict(X_test)
    roc,sens,spec=roc_auc_score(y_test,pred),recall_score(y_test,pred),recall_score(y_test,pred,pos_label=0)
    f.write(f'Recall:{sens}\tSpecificity:{spec}\tROC:{roc}\n')
    if roc>best_roc:
        best_roc=roc
        best_roc_n=n_features
        best_roc_feature_set=cols
    if sens>best_sens:
        best_sens=sens
        best_sens_n=n_features
        #best_sens_feature_set=cols
    if spec>best_spec:
        best_spec=spec
        best_spec_n=n_features
        #best_spec_feature_set=cols

    feature_imps=np.array(base_model.feature_importances_)
    non_zero_indices=np.where(feature_imps!=0)[0]
    non_zero_feature_imps=feature_imps[non_zero_indices]
    X_train=X_train[:,non_zero_indices]
    X_test=X_test[:,non_zero_indices]
    cols=cols[non_zero_indices]
    sorted_indices=np.argsort(non_zero_feature_imps)[::-1][:-int(np.floor(0.1*len(cols)))]
    X_train=X_train[:,sorted_indices]
    X_test=X_test[:,sorted_indices]
    cols=cols[sorted_indices]
    n_features=len(cols)
    iter+=1
    
f.write(f'Number of total iterations: {iter}\n')

base_cols=items['Columns']
best_indices=[]

for i in range(0,len(base_cols)):
    if base_cols[i] in best_roc_feature_set:
        best_indices.append(i)

X_best=X[:,best_indices]
specificity_scorer=make_scorer(recall_score,pos_label=0)
NPV_scorer=make_scorer(precision_score,pos_label=0)
cv=StratifiedKFold(n_splits=5,shuffle=True)
scorer={'recall':'recall','precision':'precision','specificity':specificity_scorer,'NPV':NPV_scorer,'accuracy':'accuracy','roc_auc':'roc_auc'}
cv_results=cross_validate(xgb_classifier,X_best,Y,scoring=scorer,cv=cv,return_train_score=True)
pickle.dump({'Matrix':X_best,'Phenotype':Y,'Columns':best_roc_feature_set},open(os.path.join(args.out_folder,DR_name,f'{DR_name}.pkl'),'wb'))
f.write(f'Best Specificity score and number of features: {best_spec}, {best_spec_n}.\n\
Best Sensitivity score and number of features: {best_sens},{best_sens_n}.\n\
Best ROC score and number of features: {best_roc},{best_roc_n}\n')
f.write(f'Shape of final matrix:{X_best.shape}\n')
f.write(f'CV Results: \n\n')
f.write(f"Mean Recall on Test set:{round(np.mean(cv_results['test_recall']),3)}\
    Standard Deviation : {round(np.std(cv_results['test_recall']),3)} \n")
f.write(f"Mean Specificity on Test set: {round(np.mean(cv_results['test_specificity']),3)}\
    Standard Deviation : {round(np.std(cv_results['test_specificity']),3)} \n ")
f.write(f"Mean PPV on Test set: {round(np.mean(cv_results['test_precision']),3)}\
    Standard Deviation : {round(np.std(cv_results['test_precision']),3)}. \n ")
f.write(f"Mean NPV on Test set: {round(np.mean(cv_results['test_NPV']),3)}.\
    Standard Deviation : {round(np.std(cv_results['test_NPV']),3)}. \n ")
f.write(f"Mean ROC on Test set: {round(np.mean(cv_results['test_roc_auc']),3)}\
    Standard Deviation : {round(np.std(cv_results['test_roc_auc']),3)} \n")
f.write(f"Mean Accuracy on Test set: {round(np.mean(cv_results['test_accuracy']),3)}\
    Standard Deviation: {round(np.std(cv_results['test_accuracy']),3)} \n")
for key in cv_results.keys():
    f.write(key+ ": " + str(cv_results[key])+ "\n")

f.close()

