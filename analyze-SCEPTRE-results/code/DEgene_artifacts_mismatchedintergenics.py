#6/8/26
#identifying genes that emerge as DE in a target agnostic manner; result of mixed modality intergenics

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

import warnings
warnings.simplefilter(action='ignore', category=Warning)


#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]

#cannot distinguish these since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

def get_guidecount_per_DE_gene(sample):
	SCEPTRE_results=pd.read_csv("../SCEPTRE-results/"+sample+"_SCEPTRE_results.csv")
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results=SCEPTRE_results.dropna()

	#subset to significant DE genes (FDR<0.1)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["significant"]]

	SCEPTRE_results["Direction"]=SCEPTRE_results["fold_change"].apply(lambda x: "Upregulated" if x>1 else "Downregulated")
	SCEPTRE_results["Modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])

	SCEPTRE_results=pd.DataFrame(SCEPTRE_results[["response_id","Direction","Modality"]].value_counts()).reset_index()

	#keep response IDs DE for >15 guides
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["count"]>15]

	genes_up_in_ko=SCEPTRE_results[(SCEPTRE_results["Direction"]=="Upregulated")&(SCEPTRE_results["Modality"]=="CRISPRko")]["response_id"].tolist()

	SCEPTRE_results_up_in_ko=SCEPTRE_results[SCEPTRE_results["response_id"].isin(genes_up_in_ko)]
	SCEPTRE_results_up_in_i=SCEPTRE_results[SCEPTRE_results["response_id"].isin(genes_up_in_ko)==False]

	return (SCEPTRE_results_up_in_ko,SCEPTRE_results_up_in_i)


d4_up_in_ko,d4_up_in_i=get_guidecount_per_DE_gene("d4_K562")
d7_up_in_ko,d7_up_in_i=get_guidecount_per_DE_gene("d7_K562")
d10_up_in_ko,d10_up_in_i=get_guidecount_per_DE_gene("d10_K562")
d7_A549_up_in_ko,d7_A549_up_in_i=get_guidecount_per_DE_gene("d7_A549")


def plot_guidecounts(up_in_ko_df,up_in_i_df,sample):

	(fig,ax)=plt.subplots(1,2,figsize=(7,9), layout="constrained")

	plt.suptitle("Genes with higher expression in CRISPRko-perturbed cells\n"+sample)

	up_in_ko_ko=up_in_ko_df[up_in_ko_df["Modality"]=="CRISPRko"]
	up_in_ko_i=up_in_ko_df[up_in_ko_df["Modality"]=="CRISPRi"]

	up_in_ko_ko=up_in_ko_ko.drop_duplicates(subset="response_id")
	up_in_ko_i=up_in_ko_i.drop_duplicates(subset="response_id")

	response_ids_in_both=list(set(up_in_ko_ko["response_id"].tolist()) & set(up_in_ko_i["response_id"].tolist()))
	up_in_ko_ko=up_in_ko_ko[up_in_ko_ko["response_id"].isin(response_ids_in_both)]
	up_in_ko_i=up_in_ko_i[up_in_ko_i["response_id"].isin(response_ids_in_both)]

	up_in_ko_ko=up_in_ko_ko.sort_values(by="count",ascending=False)
	up_in_ko_i = up_in_ko_i.set_index('response_id').reindex(up_in_ko_ko.set_index('response_id').index).reset_index()


	sns.barplot(data=up_in_ko_ko,y="response_id",x="count",ax=ax[0])
	ax[0].set_xlabel("# CRISPRko guides for which\ngene is upregulated")
	ax[0].set_ylabel("")

	sns.barplot(data=up_in_ko_i,y="response_id",x="count",ax=ax[1])
	ax[1].set_xlabel("# CRISPRi guides for which\ngene is downregulated")
	ax[1].set_ylabel("")

	plt.savefig("../figures/DEgene_artifacts_mismatchedintergenics/up_in_ko_"+sample+".png",bbox_inches="tight",dpi=600)

	#repeat for genes upregulated in CRISPRi perturbed cells

	(fig,ax)=plt.subplots(1,2,figsize=(7,9), layout="constrained")

	plt.suptitle("Genes with higher expression in CRISPRi-perturbed cells\n"+sample)

	up_in_i_ko=up_in_i_df[up_in_i_df["Modality"]=="CRISPRko"]
	up_in_i_i=up_in_i_df[up_in_i_df["Modality"]=="CRISPRi"]

	up_in_i_ko=up_in_i_ko.drop_duplicates(subset="response_id")
	up_in_i_i=up_in_i_i.drop_duplicates(subset="response_id")

	response_ids_in_both=list(set(up_in_i_ko["response_id"].tolist()) & set(up_in_i_i["response_id"].tolist()))
	up_in_i_ko=up_in_i_ko[up_in_i_ko["response_id"].isin(response_ids_in_both)]
	up_in_i_i=up_in_i_i[up_in_i_i["response_id"].isin(response_ids_in_both)]

	up_in_i_ko=up_in_i_ko.sort_values(by="count",ascending=False)
	up_in_i_i = up_in_i_i.set_index('response_id').reindex(up_in_i_ko.set_index('response_id').index).reset_index()


	sns.barplot(data=up_in_i_ko,y="response_id",x="count",ax=ax[0])
	ax[0].set_xlabel("# CRISPRko guides for which\ngene is downregulated")
	ax[0].set_ylabel("")

	sns.barplot(data=up_in_i_i,y="response_id",x="count",ax=ax[1])
	ax[1].set_xlabel("# CRISPRi guides for which\ngene is upregulated")
	ax[1].set_ylabel("")

	plt.savefig("../figures/DEgene_artifacts_mismatchedintergenics/up_in_i_"+sample+".png",bbox_inches="tight",dpi=600)


plot_guidecounts(d4_up_in_ko,d4_up_in_i,"Day 4 K562")
plot_guidecounts(d7_up_in_ko,d7_up_in_i,"Day 7 K562")
plot_guidecounts(d10_up_in_ko,d10_up_in_i,"Day 10 K562")
plot_guidecounts(d7_A549_up_in_ko,d7_A549_up_in_i,"Day 7 A549")
