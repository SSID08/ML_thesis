import pickle
import numpy as np 
import matplotlib.pyplot as plt
import re
import os

'''This script plots the average counts of moderate (missense variant, inframe deletion, etc.) and high impact variants\
(stop gained/lost, frameshift mutation) within cannoncial resistance associated genes in resistant and susceptible isolates\
'''

for file in os.listdir('./Updated_Data/Variant_count_Pickle_files/'):
    pickle_file=pickle.load(open(os.path.join('./Updated_Data/Variant_count_Pickle_files/',file),'rb'))
    f_name=re.sub('.pkl','',file)
    X=pickle_file['Matrix']
    pheno=pickle_file['Phenotype']
    cols=pickle_file['Column_IDs']
    my_idx = [i for i, item in enumerate(cols) if re.search('(MODERATE|HIGH)$', item)]
    out_cols=[cols[i] for i in my_idx]
    overall_freq=np.mean(X,axis=0)
    pos_pheno=np.where(pheno==1)[0]
    neg_pheno=np.where(pheno==0)[0]
    pos_matrix=X[pos_pheno,:]
    neg_matrix=X[neg_pheno,:]
    pos_freqs=np.mean(pos_matrix,axis=0)
    neg_freqs=np.mean(neg_matrix,axis=0)
    width=.35
    x=np.arange(len(out_cols))
    fig,ax=plt.subplots(figsize=(20,10))
    rects1 = ax.bar(x - width/2, pos_freqs[my_idx], width, label='Resistant')
    rects2 = ax.bar(x + width/2, neg_freqs[my_idx], width, label='Susceptible')
    ax.set_ylabel('Average of summed count',fontsize=15)
    ax.tick_params(axis='y',labelsize=12)
    ax.set_title(f'{f_name}',fontsize=20)
    ax.set_xticks(x)
    ax.set_xticklabels(out_cols,fontsize=10)
    ax.legend(fontsize=15)

    #ax.bar_label(rects1, padding=3)
    #ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.savefig(f'./Anno_Variant_Freq_plots/{f_name}.png',dpi=300)