from collections import defaultdict
import os

drug_dict=defaultdict(list)

with open('../tbdb.bed','r') as f:
    for l in f:
        chrom,pos1,pos2,gene_id,gene_name,drugs=l.split('\t')
        drugs=drugs.strip().split(',')
        for drug in drugs:
            drug_dict[drug].append([chrom,pos1,pos2,gene_id,gene_name])

print(drug_dict.keys())

for drug in drug_dict:
    print(drug)
    with open(os.path.join('../Per_Drug_bedFile',f'{drug}.bed'),'w') as f:
        for i in drug_dict[drug]:
            f.write('\t'.join(i)+ '\n')
            



    

        
