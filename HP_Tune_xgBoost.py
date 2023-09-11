import pickle
import subprocess as sp 
#from collections import defaultdict
import argparse
#from sklearn.experimental import enable_halving_search_cv
#from sklearn.ensemble import GradientBoostingClassifier
import xgboost as xgb
#from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_validate
from sklearn.metrics import precision_score,recall_score,accuracy_score,make_scorer
#from sklearn.metrics import confusion_matrix
import numpy as np
import os
import re

#genos=defaultdict(dict)

parser=argparse.ArgumentParser(description='Run hyper-parameter tuning')

parser.add_argument('--input',help="pickle input matrix and classification labels",required=True)
parser.add_argument('--out_folder',help="output folder path",required=True)
parser.add_argument('--feature_set',help="identifier for your input feature set",required=True)

args=parser.parse_args()
DR_name=re.sub('.pkl','',os.path.basename(args.input))

items=pickle.load(open(args.input,"rb"))
X=items["Matrix"]
Y=items["Phenotype"]

#print(X.shape)
#print(Y.shape)

#Prev parameters: l2_regularization=.01,max_iter=500,learning_rate=0.1,max_leaf_nodes=300,min_samples_leaf=40
xgb_classifier=xgb.XGBClassifier(subsample=.7,objective='binary:logistic',booster='gbtree',colsample_bytree=.7,learning_rate=.1,reg_lambda=1)
#param_distribution={"n_estimators":[100,300,500],"learning_rate":[0.01,0.05,0.1,0.5,1],'reg_lambda':[0.1,1,3,5,10]}
param_distribution={"n_estimators":[100,200,500,1000]}
lr_distribution={'learning_rate':[0.1,0.3,0.7,1,2]}
lambda_distribution={'reg_lambda':[0.1,0.5,1,2,5]}
#clf_GBT=GradientBoostingClassifier(n_estimators=500,learning_rae=0.05)
#X=np.array(X,dtype=object)
X_train,X_test,Y_train,y_test=train_test_split(X,Y,test_size=0.4,stratify=Y,shuffle=True)
#model=clf_GBT.fit(X_train,Y_train)
#predict=clf_GBT.predict(X_test)
#print(precision_score(y_test,y_pred=predict))
#print(recall_score(y_test,y_pred=predict))
#print(Y)
cv=StratifiedKFold(n_splits=5,shuffle=False)
tune1=GridSearchCV(xgb_classifier,param_distribution,cv=cv,scoring='roc_auc',n_jobs=-1,refit=False).fit(X_train,Y_train)
print('First tune done')
n_estimators=int(tune1.best_params_['n_estimators'])
tune2=GridSearchCV(xgb.XGBClassifier(subsample=.7,objective='binary:logistic',booster='gbtree',colsample_bytree=.7,reg_lambda=1,n_estimators=n_estimators), \
                   param_grid=lr_distribution,cv=cv,scoring='roc_auc',n_jobs=-1,refit=False).fit(X_train,Y_train)
print('Second Tune done')
lr=float(tune2.best_params_['learning_rate'])
tune3=GridSearchCV(xgb.XGBClassifier(subsample=.7,objective='binary:logistic',booster='gbtree',colsample_bytree=.7,n_estimators=n_estimators,learning_rate=lr), \
                   param_grid=lambda_distribution,cv=cv,scoring='roc_auc',n_jobs=-1,refit=True).fit(X_train,Y_train)
pred=tune3.predict(X_test)
reg_lambda=float(tune3.best_params_['reg_lambda'])
specificity_scorer=make_scorer(recall_score,pos_label=0)
NPV_scorer=make_scorer(precision_score,pos_label=0)
spec=recall_score(y_test,pred,pos_label=0)
NPV=precision_score(y_test,pred,pos_label=0)
recall=recall_score(y_test,pred)
precision=precision_score(y_test,pred)
acc_score=accuracy_score(y_test,pred)
final_clf=xgb.XGBClassifier(objective='binary:logistic',booster='gbtree',reg_lambda=reg_lambda,n_estimators=n_estimators,learning_rate=lr)
scorer={'recall':'recall','precision':'precision','specificity':specificity_scorer,'NPV':NPV_scorer,'accuracy':'accuracy','roc_auc':'roc_auc'}
cv_results=cross_validate(final_clf,X,Y,scoring=scorer,cv=StratifiedKFold(n_splits=5,shuffle=True),return_train_score=True)
#print(str(param_distribution))
print(DR_name)
print(f'Best CV score: {str(tune2.best_score_)}; Best N_estimators: {str(n_estimators)}; Best_Learning_Rate : {str(lr)}')
print(f'Recall score on validation set: {recall}')
print(f'Specificity score on validation set: {spec}')
print(f'NPV on validation set: {NPV}')
print(f'PPV score on validation set : {precision}')
print(f'Accuracy score on validation set: {acc_score}')

#pickle.dump({"model":tune,"Recall_test":str(recall),"Precision_test":str(precision)}, open(os.path.join(args.out_folder,f"xgBoost_{DR_name}_{args.feature_set}_"+"tune_model.pkl"),"wb"))
with open(os.path.join(args.out_folder,DR_name,f"xgBoost_{DR_name}_{args.feature_set}"+".txt"),'w') as f:
    f.write(f'Best CV score: {str(tune3.best_score_)}; Best N_estimators: {str(n_estimators)}; Best_Learning_Rate : {str(lr)}; Best reg_lamba={str(reg_lambda)}\n')
    f.write(f'Recall score on validation set: {recall}\n Specificity score on validation set {spec}\n')
    f.write(f'NPV on validation set: {NPV}\nPPV score on validation set : {precision}\n')
    f.write(f'Accuracy score on test set: {acc_score}\n')
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