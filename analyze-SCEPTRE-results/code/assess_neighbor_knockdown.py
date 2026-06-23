#6/1/2026
#do we see knockdown for genes whose MANE select TSS is within 1kb of another TSS? 

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import warnings
warnings.simplefilter(action='ignore', category=Warning)

# -- Reading in, transforming data --

#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")

guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]
guide_to_shortid=guideref.set_index("guide_id_long").to_dict()["guide_id_short"]

#cannot distinguish these at Day 4,7,10 since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]


def processdata_nonday14(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna()
	SCEPTRE_results["grna_id_short"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_shortid[x])
	SCEPTRE_results=SCEPTRE_results.sort_values(by="grna_id_short")
	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]

	SCEPTRE_results_CRISPRko_target=SCEPTRE_results_CRISPRko[SCEPTRE_results_CRISPRko["response_id"]==SCEPTRE_results_CRISPRko["grna_target"]].reset_index()
	guides_with_target_detected=SCEPTRE_results_CRISPRko_target["grna_id"].unique() 
	SCEPTRE_results_CRISPRko=SCEPTRE_results_CRISPRko[SCEPTRE_results_CRISPRko["grna_id"].isin(guides_with_target_detected)].reset_index(drop=True)

	SCEPTRE_results_CRISPRi_target=SCEPTRE_results_CRISPRi[SCEPTRE_results_CRISPRi["response_id"]==SCEPTRE_results_CRISPRi["grna_target"]].reset_index()
	guides_with_target_detected=SCEPTRE_results_CRISPRi_target["grna_id"].unique() 
	SCEPTRE_results_CRISPRi=SCEPTRE_results_CRISPRi[SCEPTRE_results_CRISPRi["grna_id"].isin(guides_with_target_detected)].reset_index(drop=True)

	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi,SCEPTRE_results_CRISPRko_target,SCEPTRE_results_CRISPRi_target)

day7_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRi,_,_=processdata_nonday14("../SCEPTRE-results/d7_K562_SCEPTRE_results.csv")

#identify instances of neighbor KD
targets_with_neighboring_genes=pd.read_csv("../../reference/genes_with_MANEselectTSS_within1kb.csv")
targets_with_neighboring_genes=targets_with_neighboring_genes[targets_with_neighboring_genes["Selected Gene name"].isin(day7_SCEPTRE_results_CRISPRi["grna_target"].unique().tolist())]
target_to_neighboring_tss=targets_with_neighboring_genes.set_index("Selected Gene name").to_dict()["Overlapping TSS"]
target_to_neighboring_genename=targets_with_neighboring_genes.set_index("Selected Gene name").to_dict()["Overlapping Gene name"]

day7_SCEPTRE_results_CRISPRi=day7_SCEPTRE_results_CRISPRi[day7_SCEPTRE_results_CRISPRi["grna_target"].isin(targets_with_neighboring_genes["Selected Gene name"].tolist())]
day7_SCEPTRE_results_CRISPRi["Neighboring gene"]=day7_SCEPTRE_results_CRISPRi["grna_target"].apply(lambda x:target_to_neighboring_genename[x])
day7_SCEPTRE_results_CRISPRi["Neighboring TSS position"]=day7_SCEPTRE_results_CRISPRi["grna_target"].apply(lambda x:target_to_neighboring_tss[x])
day7_SCEPTRE_results_CRISPRi["guide sequence"]=day7_SCEPTRE_results_CRISPRi["grna_id"].apply(lambda x:x.split("_")[1])

#get guide "cut" positions: 17nt along the 20nt guide sequence
CRISPick_results_CRISPRi=pd.read_table("../../reference/CRISPRi_Cas9_CROPseqbenchmark-sgrna-designs.txt")
guidesequence_to_position=CRISPick_results_CRISPRi.set_index("sgRNA Sequence").to_dict()["sgRNA 'Cut' Position"]
day7_SCEPTRE_results_CRISPRi["Guide position"]=day7_SCEPTRE_results_CRISPRi["guide sequence"].apply(lambda x: guidesequence_to_position[x])
day7_SCEPTRE_results_CRISPRi["Distance to neighboring TSS"]=day7_SCEPTRE_results_CRISPRi.apply(lambda x:int(abs(x["Guide position"]-x["Neighboring TSS position"])),axis=1)

day7_SCEPTRE_results_CRISPRi["Guide index"]=day7_SCEPTRE_results_CRISPRi["grna_id_short"].apply(lambda x:x.split("_")[1])

day7_SCEPTRE_results_CRISPRi=day7_SCEPTRE_results_CRISPRi.sort_values(by="Distance to neighboring TSS")

#violin plot in LFC of all transcripts, blue dot indicates relative position of target. Day 7
(fig,ax)=plt.subplots(2,len(target_to_neighboring_genename.keys()),
	figsize=(len(target_to_neighboring_genename.keys())*2,4), 
	layout="constrained",
	height_ratios=[6, 1],
	width_ratios=[4, 1,4,2,4],
	sharey="row")
plt.suptitle("Knockdown of genes with TSS within 1kb of target (per guide)\nDay7 CRISPRi")

axes_idx=0

for target in day7_SCEPTRE_results_CRISPRi["grna_target"].unique().tolist():
	SCEPTRE_subset_to_target=day7_SCEPTRE_results_CRISPRi[day7_SCEPTRE_results_CRISPRi["grna_target"]==target]
	neighboring_gene=target_to_neighboring_genename[target]
	neighboring_gene_data=SCEPTRE_subset_to_target[SCEPTRE_subset_to_target["response_id"]==neighboring_gene]

	sns.violinplot(data=SCEPTRE_subset_to_target,x="Guide index",y="log(foldchange)",inner=None,color="lightgrey",ax=ax[0,axes_idx])
	ax[0,axes_idx].margins(x=0.01)
	ax[0,axes_idx].set_xlabel("")
	ax[0,axes_idx].set_ylabel("")
	ax[0,axes_idx].set_title("Target: "+target,fontsize=12)
	ax[0,axes_idx].scatter(data=neighboring_gene_data[neighboring_gene_data["significant"]==False],x="Guide index",y="log(foldchange)",color="blue",s=125)
	ax[0,axes_idx].scatter(data=neighboring_gene_data[neighboring_gene_data["significant"]],x="Guide index",y="log(foldchange)",color="blue",s=225,marker='*')

	ax[1,axes_idx].set_title("bp to "+neighboring_gene+":")

	sns.heatmap(data=SCEPTRE_subset_to_target[["Distance to neighboring TSS"]].drop_duplicates(keep="first").T,
		annot=True, fmt="",
		cmap="Blues", vmin=0, vmax=10000,cbar=False,
		ax=ax[1,axes_idx],
		alpha=0)
	ax[1,axes_idx].set_xticks([])
	ax[1,axes_idx].set_yticks([])
	ax[1,axes_idx].set(xlabel=None, ylabel=None)

	axes_idx=axes_idx+1

ax[0,0].set_ylabel("log(foldchange)")

plt.savefig("../figures/day7_CRISPRi_neighbordownreg_violinplots.png",bbox_inches="tight",dpi=600)

#repeat for CRISPRko

day7_SCEPTRE_results_CRISPRko["Guide index"]=day7_SCEPTRE_results_CRISPRko["grna_id_short"].apply(lambda x:x.split("_")[1])

(fig,ax)=plt.subplots(1,len(target_to_neighboring_genename.keys()),
	figsize=(len(target_to_neighboring_genename.keys())*2,3), 
	layout="constrained",
	width_ratios=[4, 2,4,3,4],
	sharey=True)

plt.suptitle("Knockdown of genes with TSS within 1kb of target (per guide)\nDay7 CRISPRko")
axes_idx=0

#this is intentionally the same order as CRISPRi
for target in day7_SCEPTRE_results_CRISPRi["grna_target"].unique().tolist():
	SCEPTRE_subset_to_target=day7_SCEPTRE_results_CRISPRko[day7_SCEPTRE_results_CRISPRko["grna_target"]==target]
	neighboring_gene=target_to_neighboring_genename[target]
	neighboring_gene_data=SCEPTRE_subset_to_target[SCEPTRE_subset_to_target["response_id"]==neighboring_gene]

	sns.violinplot(data=SCEPTRE_subset_to_target,x="Guide index",y="log(foldchange)",inner=None,color="lightgrey",ax=ax[axes_idx])
	ax[axes_idx].margins(x=0.01)
	ax[axes_idx].set_xlabel("")
	ax[axes_idx].set_ylabel("")
	ax[axes_idx].set_title("Target: "+target+"\nNeighbor: "+neighboring_gene,fontsize=11)
	ax[axes_idx].scatter(data=neighboring_gene_data[neighboring_gene_data["significant"]==False],x="Guide index",y="log(foldchange)",color="blue",s=125)
	ax[axes_idx].scatter(data=neighboring_gene_data[neighboring_gene_data["significant"]],x="Guide index",y="log(foldchange)",color="blue",s=225,marker='*')

	axes_idx=axes_idx+1
ax[0].set_ylabel("log(foldchange)")



plt.savefig("../figures/day7_CRISPRko_neighbordownreg_violinplots.png",bbox_inches="tight",dpi=600)

