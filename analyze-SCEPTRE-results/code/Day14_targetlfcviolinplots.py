#5/29/26
#identifying log fold change of target knockdown from SCEPTRE results relative to all transcripts

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

import warnings
warnings.simplefilter(action='ignore', category=Warning)


#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_shortid=guideref.set_index("guide_id_long").to_dict()["guide_id_short"]

def process_SCEPTRE_results(modality):
	SCEPTRE_results=pd.read_csv("../SCEPTRE-results/d14_K562_"+modality+"_SCEPTRE_results.csv")
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna()
	SCEPTRE_results["grna_id_short"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_shortid[x])
	SCEPTRE_results=SCEPTRE_results.sort_values(by="grna_id_short")

	SCEPTRE_results_target=SCEPTRE_results[SCEPTRE_results["response_id"]==SCEPTRE_results["grna_target"]].reset_index()
	guides_with_target_detected=SCEPTRE_results["grna_id"].unique() 
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_with_target_detected)].reset_index(drop=True)


	return SCEPTRE_results,SCEPTRE_results_target

SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRko_target=process_SCEPTRE_results("CRISPRko")
SCEPTRE_results_CRISPRi,SCEPTRE_results_CRISPRi_target=process_SCEPTRE_results("CRISPRi")


(fig,ax)=plt.subplots(2,1,figsize=(11,6), layout="constrained")

sns.violinplot(data=SCEPTRE_results_CRISPRko,x="grna_id_short",y="log(foldchange)",inner=None,color="lightgrey",ax=ax[0])
ax[0].tick_params(axis='x', rotation=90)
ax[0].margins(x=0.01)
ax[0].set_xlabel("")
ax[0].set_ylabel("log(fold change + 1)")
ax[0].set_title("Day 14 CRISPRko")
SCEPTRE_results_CRISPRko_target_notsig=SCEPTRE_results_CRISPRko_target[SCEPTRE_results_CRISPRko_target["significant"]==False]
ax[0].scatter(SCEPTRE_results_CRISPRko_target_notsig.index,SCEPTRE_results_CRISPRko_target_notsig["log(foldchange)"],label="Target transcript",color="blue",s=30)
SCEPTRE_results_CRISPRko_target_sig=SCEPTRE_results_CRISPRko_target[SCEPTRE_results_CRISPRko_target["significant"]]
ax[0].scatter(SCEPTRE_results_CRISPRko_target_sig.index,SCEPTRE_results_CRISPRko_target_sig["log(foldchange)"],label="Target transcript",color="blue",s=95,marker='*')

sns.violinplot(data=SCEPTRE_results_CRISPRi,x="grna_id_short",y="log(foldchange)",inner=None,color="lightgrey",ax=ax[1])
ax[1].tick_params(axis='x', rotation=90)
ax[1].margins(x=0.01)
ax[1].set_xlabel("")
ax[1].set_ylabel("log(fold change + 1)")
ax[1].set_title("Day 14 CRISPRi")
SCEPTRE_results_CRISPRi_target_notsig=SCEPTRE_results_CRISPRi_target[SCEPTRE_results_CRISPRi_target["significant"]==False]
ax[1].scatter(SCEPTRE_results_CRISPRi_target_notsig.index,SCEPTRE_results_CRISPRi_target_notsig["log(foldchange)"],label="Target transcript",color="blue",s=30)
SCEPTRE_results_CRISPRi_target_sig=SCEPTRE_results_CRISPRi_target[SCEPTRE_results_CRISPRi_target["significant"]]
ax[1].scatter(SCEPTRE_results_CRISPRi_target_sig.index,SCEPTRE_results_CRISPRi_target_sig["log(foldchange)"],label="Target transcript",color="blue",s=95,marker='*')

fig.suptitle("Perturbation effect on all qc-passing transcripts, target in blue")

plt.savefig("../figures/K562_Day14_Target_Downregulation_violinplots.png",bbox_inches="tight",dpi=600)


