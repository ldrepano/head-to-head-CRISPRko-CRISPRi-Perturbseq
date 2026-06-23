#5/29/2026
#identifying promiscuous guides; ~=those poorly correlated with others targeting the same gene

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import warnings
warnings.simplefilter(action='ignore', category=Warning)
from adjustText import adjust_text

#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]
guide_to_shortid=guideref.set_index("guide_id_long").to_dict()["guide_id_short"]


#cannot distinguish these at Day 4,7,10 since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

def processdata(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results=SCEPTRE_results.dropna()

	#creates a column in the supplied df called signed_logpval
	#represents each guide/transcript pair such that magnitude indicates significance of perturbation effect relative to negative controls
	#... and sign (negative vs positive) indicates upregulation vs downregulation
	SCEPTRE_results["effect_sign"]=SCEPTRE_results["fold_change"].apply(lambda x: -1 if x<1 else 1)
	SCEPTRE_results["-logpval"]=-np.log(SCEPTRE_results["p_value"])
	SCEPTRE_results["signed_logpval"]=SCEPTRE_results["effect_sign"]*SCEPTRE_results["-logpval"]
	SCEPTRE_results=SCEPTRE_results.drop(["effect_sign","-logpval"],axis=1)

	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]

	#transform to one row per response gene, one column per guide
	SCEPTRE_results_CRISPRko=pd.pivot_table(SCEPTRE_results_CRISPRko[["response_id","grna_id","signed_logpval"]], values="signed_logpval", index="response_id", columns=["grna_id"])
	SCEPTRE_results_CRISPRi=pd.pivot_table(SCEPTRE_results_CRISPRi[["response_id","grna_id","signed_logpval"]], values="signed_logpval", index="response_id", columns=["grna_id"])


	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi)

day4_SCEPTRE_results_CRISPRko,day4_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d4_K562_SCEPTRE_results.csv")
day7_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d7_K562_SCEPTRE_results.csv")
day10_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d10_K562_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRko,_=processdata("../SCEPTRE-results/d14_K562_CRISPRko_SCEPTRE_results.csv")
_,day14_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d14_K562_CRISPRi_SCEPTRE_results.csv")



def get_correlations(timepoint1,timepoint2):

	#identify correlation of each guide with itself across replicates
	responsegenes_in_both=list(set(timepoint1.index.tolist()).intersection(set(timepoint2.index.tolist())))
	timepoint1=timepoint1.loc[responsegenes_in_both]
	timepoint1=timepoint1.sort_index()
	timepoint2=timepoint2.loc[responsegenes_in_both]
	timepoint2=timepoint2.sort_index()


	crosstimecorr=timepoint1.corrwith(timepoint2,numeric_only=True)
	crosstimecorr=crosstimecorr.dropna()
	guide_to_crosstime_correlation=pd.DataFrame(crosstimecorr).to_dict()[0]


	#get correlation between guides targeting the same gene (timepoint1)
	timepoint1 = timepoint1.corr(numeric_only=True)
	timepoint1_long=pd.melt(timepoint1,ignore_index=False)
	timepoint1_long["other_grna_id"]=timepoint1_long.index
	timepoint1_long=timepoint1_long[timepoint1_long["grna_id"]!=timepoint1_long["other_grna_id"]]
	timepoint1_long["same_target_gene"]=timepoint1_long.apply(lambda x: x["grna_id"].split("_")[0]==x["other_grna_id"].split("_")[0],axis=1)
	timepoint1_samegene=timepoint1_long[timepoint1_long["same_target_gene"]]

	#repeat for rep2
	timepoint2 = timepoint2.corr(numeric_only=True)
	timepoint2_long=pd.melt(timepoint2,ignore_index=False)
	timepoint2_long["other_grna_id"]=timepoint2_long.index
	timepoint2_long=timepoint2_long[timepoint2_long["grna_id"]!=timepoint2_long["other_grna_id"]]
	timepoint2_long["same_target_gene"]=timepoint2_long.apply(lambda x: x["grna_id"].split("_")[0]==x["other_grna_id"].split("_")[0],axis=1)
	timepoint2_samegene=timepoint2_long[timepoint2_long["same_target_gene"]]

	#average correlation between guides targeting same gene across time points
	samegene=pd.concat([timepoint1_samegene,timepoint2_samegene])
	samegene_mediancorrbyguide=samegene.groupby("other_grna_id").agg(median_corr=("value",'median')).reset_index()

	samegene_mediancorrbyguide=samegene_mediancorrbyguide[samegene_mediancorrbyguide["other_grna_id"].isin(guide_to_crosstime_correlation.keys())]
	samegene_mediancorrbyguide["guide_crosstime_corr"]=samegene_mediancorrbyguide["other_grna_id"].apply(lambda x:guide_to_crosstime_correlation[x])

	return samegene_mediancorrbyguide



d4_v_7_CRISPRko_correlations=get_correlations(day4_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRko)
d4_v_7_CRISPRi_correlations=get_correlations(day4_SCEPTRE_results_CRISPRi,day7_SCEPTRE_results_CRISPRi)


d7_v_10_CRISPRko_correlations=get_correlations(day7_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRko)
d7_v_10_CRISPRi_correlations=get_correlations(day7_SCEPTRE_results_CRISPRi,day10_SCEPTRE_results_CRISPRi)

d14_v_10_CRISPRko_correlations=get_correlations(day10_SCEPTRE_results_CRISPRko,day14_SCEPTRE_results_CRISPRko)
d14_v_10_CRISPRi_correlations=get_correlations(day10_SCEPTRE_results_CRISPRi,day14_SCEPTRE_results_CRISPRi)


def make_correlation_scatterplot(df,name):
	plt.figure()
	plt.scatter(data=df,x="median_corr",y="guide_crosstime_corr",c="black")

	#label guides whose replicate correlation is at least 3x its correlation to other guides targeting the same gene
	potential_promiscuous_guides=df[(df["median_corr"]*3) < (df["guide_crosstime_corr"])]
	promiscuous_guide_labels=[guide_to_shortid[guide] for guide in potential_promiscuous_guides["other_grna_id"].tolist()]
	promiscuous_guide_targetcorr=potential_promiscuous_guides["median_corr"].tolist()
	texts=[]
	promiscuous_guide_repcorr=potential_promiscuous_guides["guide_crosstime_corr"].tolist()
	for i, txt in enumerate(promiscuous_guide_labels):
		texts.append(plt.text(promiscuous_guide_targetcorr[i], promiscuous_guide_repcorr[i], txt,c="blue"))
		
	adjust_text(texts, arrowprops=dict(arrowstyle='->', color='blue'),force_text=(0.2,0.4))

	plt.xlabel("Median correlation of guide with others targeting the same gene")
	plt.ylabel("Guide correlation across time points\n(separate transduction)")
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.title(name)
	plt.axline([0, 0], slope=1, color='black', linestyle='--')
	plt.savefig("../figures/"+name+"_vs_withingene_correlation_K562.png",bbox_inches="tight",dpi=600)


make_correlation_scatterplot(d4_v_7_CRISPRko_correlations,"Day 4 vs 7 CRISPRko")
make_correlation_scatterplot(d4_v_7_CRISPRi_correlations,"Day 4 vs 7 CRISPRi")
make_correlation_scatterplot(d7_v_10_CRISPRko_correlations,"Day 7 vs 10 CRISPRko")
make_correlation_scatterplot(d7_v_10_CRISPRi_correlations,"Day 7 vs 10 CRISPRi")
make_correlation_scatterplot(d14_v_10_CRISPRko_correlations,"Day 10 vs 14 CRISPRko")
make_correlation_scatterplot(d14_v_10_CRISPRi_correlations,"Day 10 vs 14 CRISPRi")



