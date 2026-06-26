#6/8/26
#identifying transcriptional effect size of perturbing each gene in library 

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
plt.rc('pdf', fonttype=42)


import warnings
warnings.simplefilter(action='ignore', category=Warning)


#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]

#cannot distinguish these since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

def get_n_de_genes_per_guide(sample):
	SCEPTRE_results=pd.read_csv("../SCEPTRE-results/"+sample+"_SCEPTRE_results.csv")
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results=SCEPTRE_results.dropna()

	#subset to significant DE genes (FDR<0.1)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["significant"]]

	SCEPTRE_results["Target"]=SCEPTRE_results["grna_id"].apply(lambda x:x.split("_")[0])

	SCEPTRE_results["Essential positive control?"]=SCEPTRE_results["Target"].isin(["RPL9","DCAF13"])


	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]

	SCEPTRE_results_CRISPRko=pd.DataFrame(SCEPTRE_results_CRISPRko[["grna_id","Target","Essential positive control?"]].value_counts()).reset_index()
	SCEPTRE_results_CRISPRi=pd.DataFrame(SCEPTRE_results_CRISPRi[["grna_id","Target","Essential positive control?"]].value_counts()).reset_index()
	
	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi)


d4_CRISPRko,d4_CRISPRi=get_n_de_genes_per_guide("d4_K562")
d7_CRISPRko,d7_CRISPRi=get_n_de_genes_per_guide("d7_K562")


(fig,ax)=plt.subplots(1,2,figsize=(8,3), layout="constrained")

mean_order = d4_CRISPRko.groupby('Target')['count'].mean().sort_values(ascending=False)
d4_CRISPRko = d4_CRISPRko.set_index('Target').loc[mean_order.index].reset_index()

sns.barplot(data=d4_CRISPRko,x="Target",y="count",hue="Essential positive control?",ax=ax[0],palette=["lightgrey","blue"],errorbar=None)
sns.stripplot(data=d4_CRISPRko,x="Target", y="count",ax=ax[0],c="black",jitter=0,size=3)
ax[0].set_ylabel("# affected genes")
ax[0].tick_params(axis='x', labelrotation=90)
ax[0].set_title("CRISPRko Day 4")

mean_order = d4_CRISPRi.groupby('Target')['count'].mean().sort_values(ascending=False)
d4_CRISPRi = d4_CRISPRi.set_index('Target').loc[mean_order.index].reset_index()

sns.barplot(data=d4_CRISPRi,x="Target",y="count",hue="Essential positive control?",ax=ax[1],palette=["lightgrey","blue"],errorbar=None)
sns.stripplot(data=d4_CRISPRi,x="Target", y="count",ax=ax[1],c="black",jitter=0,size=3)

ax[1].set_ylabel("# affected genes")
ax[1].tick_params(axis='x', labelrotation=90)
ax[1].set_title("CRISPRi Day 4")

plt.savefig("../figures/n_DE_genes_per_target_d4.pdf",bbox_inches="tight",dpi=600)

(fig,ax)=plt.subplots(1,2,figsize=(8,3), layout="constrained")

mean_order = d7_CRISPRko.groupby('Target')['count'].mean().sort_values(ascending=False)
d7_CRISPRko = d7_CRISPRko.set_index('Target').loc[mean_order.index].reset_index()

sns.barplot(data=d7_CRISPRko,x="Target",y="count",hue="Essential positive control?",ax=ax[0],palette=["lightgrey","blue"],errorbar=None)
sns.stripplot(data=d7_CRISPRko,x="Target", y="count",ax=ax[0],c="black",jitter=0,size=3)
ax[0].set_ylabel("# affected genes")
ax[0].tick_params(axis='x', labelrotation=90)
ax[0].set_title("CRISPRko Day 7")

mean_order = d7_CRISPRi.groupby('Target')['count'].mean().sort_values(ascending=False)
d7_CRISPRi = d7_CRISPRi.set_index('Target').loc[mean_order.index].reset_index()

sns.barplot(data=d7_CRISPRi,x="Target",y="count",hue="Essential positive control?",ax=ax[1],palette=["lightgrey","blue"],errorbar=None)
sns.stripplot(data=d7_CRISPRi,x="Target", y="count",ax=ax[1],c="black",jitter=0,size=3)
ax[1].set_ylabel("# affected genes")
ax[1].tick_params(axis='x', labelrotation=90)
ax[1].set_title("CRISPRi Day 7")

plt.savefig("../figures/n_DE_genes_per_target_d7.pdf",bbox_inches="tight",dpi=600)

