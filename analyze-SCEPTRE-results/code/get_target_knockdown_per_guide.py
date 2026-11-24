#11/17/26
#saving target DE results for all samples (Supplementary Data 3)

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=Warning)
from adjustText import adjust_text
plt.rc('pdf', fonttype=42)


color_dictionary={"CRISPRko": 'dodgerblue', "CRISPRi": 'limegreen'}

#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]
guide_to_shortid=guideref.set_index("guide_id_long").to_dict()["guide_id_short"]


#cannot distinguish these at Day 4,7,10 since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

def get_target_kd_per_guide(sample):
    print(sample)
    #get SCEPTRE results
    SCEPTRE_results=pd.read_csv("../SCEPTRE-results/"+sample+"_SCEPTRE_results.csv")
    SCEPTRE_results=SCEPTRE_results.dropna()
    SCEPTRE_results['sample'] = sample

    if 'd14' not in sample:
        SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
        SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
    else:
        SCEPTRE_results['modality']=sample.split('_')[2]
        SCEPTRE_results['sample']=sample.split('_')[0]+'_'+sample.split('_')[1]
        
    
    #get target lfc and pvalue
    SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log10(x+0.1))
    SCEPTRE_target_results = SCEPTRE_results[SCEPTRE_results['response_id']==SCEPTRE_results['grna_target']][['sample',
                                                                                                              'modality',
                                                                                                              'grna_id',
                                                                                                              'grna_target',
                                                                                                              'log(foldchange)',
                                                                                                              'p_value',
                                                                                                              'significant']].drop_duplicates()
    return SCEPTRE_target_results

samplenames = ["d4_K562","d7_K562","d10_K562","d14_K562_CRISPRko","d14_K562_CRISPRi","d7_A549"]
target_kd_df = pd.concat([get_target_kd_per_guide(samplename) for samplename in samplenames]).reset_index(drop=True)
target_kd_df["guide_id"]=target_kd_df["grna_id"].apply(lambda x: guide_to_shortid[x])

target_kd_df.to_csv('../SCEPTRE-results/target_SCEPTRE_results.csv',index=False)