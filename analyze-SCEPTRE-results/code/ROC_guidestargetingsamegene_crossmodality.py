#6/8/2026
#Day 14 K562 examine guide agreement across CRISPRko and CRISPRi, with and without promiscuous CRISPRi guides removed

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import warnings
warnings.simplefilter(action='ignore', category=Warning)
plt.rc('pdf', fonttype=42)



#guides identified as promiscuous from cross-guide vs cross-sample correlation analysis 
promiscuous_CRISPRi_labels=["STK24_3","SH3KBP1_4","PYGB_2","IFT25_3","AKT2_2", "TMSB10_1","PDHA1_3"]
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
short_to_long_id=guideref.set_index("guide_id_short").to_dict()["guide_id_long"]
promiscuous_CRISPRi_guides=[short_to_long_id[guide] for guide in promiscuous_CRISPRi_labels]

def processdata(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna()
	return SCEPTRE_results


day14_SCEPTRE_results_CRISPRko=processdata("../SCEPTRE-results/d14_K562_CRISPRko_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d14_K562_CRISPRi_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRko["modality"]="CRISPRko"
day14_SCEPTRE_results_CRISPRi["modality"]="CRISPRi"
day14_SCEPTRE_results=pd.concat([day14_SCEPTRE_results_CRISPRko,day14_SCEPTRE_results_CRISPRi])


def get_signed_negative_log_pval(SCEPTRE_results_df):
	#creates a column in the supplied df called signed_logpval
	#represents each guide/transcript pair such that magnitude indicates significance of perturbation effect relative to negative controls
	#... and sign (negative vs positive) indicates upregulation vs downregulation
	SCEPTRE_results_df["effect_sign"]=SCEPTRE_results_df["fold_change"].apply(lambda x: -1 if x<1 else 1)
	SCEPTRE_results_df["-logpval"]=-np.log(SCEPTRE_results_df["p_value"])
	SCEPTRE_results_df["signed_logpval"]=SCEPTRE_results_df["effect_sign"]*SCEPTRE_results_df["-logpval"]
	SCEPTRE_results_df=SCEPTRE_results_df.drop(["effect_sign","-logpval"],axis=1)
	return SCEPTRE_results_df

#get dataframe of true positive rates, false positive rates for all signed p-value cutoffs to distinguish guides targeting same gene
def get_tpr_fpr_df(SCEPTRE_results,guides_included,remove_promiscuous_guides=False):
	SCEPTRE_results=get_signed_negative_log_pval(SCEPTRE_results)
	signedlogpvals=SCEPTRE_results[["response_id","grna_id","signed_logpval","modality"]]
	signedlogpvals=pd.pivot_table(signedlogpvals, values="signed_logpval", index="response_id", columns=["grna_id","modality"]).reset_index()
	#get pearson correlation in signed pvals attached to response gene for each guide 
	correlations = signedlogpvals.corr(numeric_only=True)

	#ROC-AUC of guides targeting the same vs different gene 
	#remove self-correlation by excluding the first guide from the other set 
	correlations_long=pd.melt(correlations,ignore_index=False)
	correlations_long=correlations_long.rename(columns={"grna_id":"guide2","modality":"modality2"})
	correlations_long=correlations_long.reset_index()
	#remove those of the same modality by requiring guide1 to be CRISPRko, guide2 to be CRISPRi 
	correlations_long=correlations_long[correlations_long["modality"]=="CRISPRko"]
	correlations_long=correlations_long[correlations_long["modality2"]=="CRISPRi"]

	if remove_promiscuous_guides:
		correlations_long=correlations_long[correlations_long["guide2"].isin(promiscuous_CRISPRi_guides)==False]


	correlations_long["same_target_gene"]=correlations_long.apply(lambda x: x["grna_id"].split("_")[0]==x["guide2"].split("_")[0],axis=1)
	fpr, tpr, threshold = roc_curve(correlations_long['same_target_gene'], correlations_long['value'])
	tpr_fpr_df = pd.DataFrame({'True Positive Rate': tpr, 'False Positive Rate': fpr, 'threshold': threshold})
	roc_auc = auc(fpr, tpr)
	tpr_fpr_df["Guides included"]=guides_included+", AUC="+str(round(roc_auc,2))
	return tpr_fpr_df

tpr_fpr_df_d14=get_tpr_fpr_df(day14_SCEPTRE_results,guides_included="All")
tpr_fpr_df_d14_nopromiscuous=get_tpr_fpr_df(day14_SCEPTRE_results,guides_included="Nonpromiscuous",remove_promiscuous_guides=True)

tpr_fpr_df=pd.concat([tpr_fpr_df_d14_nopromiscuous,tpr_fpr_df_d14])

plt.subplots(figsize=(4, 4))
sns.lineplot(data=tpr_fpr_df, x='False Positive Rate',
             y='True Positive Rate',style="Guides included",errorbar=None,color="black")
sns.despine()
plt.title("Cross-modality agreement in perturbation effects\nDay 14 K562")
plt.savefig("../figures/crossmodality_ROC_signedpvalcorrelation_d14K562.pdf",bbox_inches="tight",dpi=600)


