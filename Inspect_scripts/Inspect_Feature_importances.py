import numpy as np 
import os
import pickle

model=pickle.load(open('./Final_Results/amikacin/Models/XGB_Boost_amikacin_All_DR.pkl',"rb"))
features=pickle.load(open('./Files/Transformed_dataframes_All_DR/amikacin.pkl','rb'))['Columns']
model=model['Model']
feature_imp=model.feature_importances_
sorted_indices=np.argsort(feature_imp)[::-1][0:20]
top_20_features=feature_imp[sorted_indices]
top_20_feature_names=features[sorted_indices]
print(np.sort(model.feature_importances_)[::-1][0:10])