#5/29/2026

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import warnings
warnings.simplefilter(action='ignore', category=Warning)
plt.rc('pdf', fonttype=42)


#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]


#cannot distinguish these at Day 4,7,10 since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

def processdata_nonday14(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna()
	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]
	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi)

day4_SCEPTRE_results_CRISPRko,day4_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d4_K562_SCEPTRE_results.csv")
day7_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d7_K562_SCEPTRE_results.csv")
day10_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d10_K562_SCEPTRE_results.csv")


day14_SCEPTRE_results_CRISPRko=pd.read_csv("../SCEPTRE-results/d14_K562_CRISPRko_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRko["log(foldchange)"]=day14_SCEPTRE_results_CRISPRko["fold_change"].apply(lambda x: np.log(x+0.1))
day14_SCEPTRE_results_CRISPRko=day14_SCEPTRE_results_CRISPRko.dropna()

day14_SCEPTRE_results_CRISPRi=pd.read_csv("../SCEPTRE-results/d14_K562_CRISPRi_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRi["log(foldchange)"]=day14_SCEPTRE_results_CRISPRi["fold_change"].apply(lambda x: np.log(x+0.1))
day14_SCEPTRE_results_CRISPRi=day14_SCEPTRE_results_CRISPRi.dropna()

def get_signed_negative_log10_pval(SCEPTRE_results_df):
	#creates a column in the supplied df called signed_logpval
	#represents each guide/transcript pair such that magnitude indicates significance of perturbation effect relative to negative controls
	#... and sign (negative vs positive) indicates upregulation vs downregulation
	SCEPTRE_results_df["effect_sign"]=SCEPTRE_results_df["fold_change"].apply(lambda x: -1 if x<1 else 1)
	SCEPTRE_results_df["-logpval"]=-np.log(SCEPTRE_results_df["p_value"])
	SCEPTRE_results_df["signed_logpval"]=SCEPTRE_results_df["effect_sign"]*SCEPTRE_results_df["-logpval"]
	SCEPTRE_results_df=SCEPTRE_results_df.drop(["effect_sign","-logpval"],axis=1)
	return SCEPTRE_results_df

#get dataframe of true positive rates, false positive rates for all signed p-value cutoffs to distinguish guides targeting same gene
def get_tpr_fpr_df(SCEPTRE_results,day):
	SCEPTRE_results=get_signed_negative_log10_pval(SCEPTRE_results)
	signedlogpvals=SCEPTRE_results[["response_id","grna_id","signed_logpval"]]
	signedlogpvals=pd.pivot_table(signedlogpvals, values="signed_logpval", index="response_id", columns=["grna_id"]).reset_index()
	#get pearson correlation in signed pvals attached to response gene for each guide 
	correlations = signedlogpvals.corr(numeric_only=True)

	#ROC-AUC of guides targeting the same vs different gene 
	#remove self-correlation by excluding the first guide from the other set 
	correlations_long=pd.melt(correlations,ignore_index=False)
	correlations_long["other_grna_id"]=correlations_long.index
	correlations_long=correlations_long[correlations_long["grna_id"]!=correlations_long["other_grna_id"]]
	correlations_long["same_target_gene"]=correlations_long.apply(lambda x: x["grna_id"].split("_")[0]==x["other_grna_id"].split("_")[0],axis=1)
	fpr, tpr, threshold = roc_curve(correlations_long['same_target_gene'], correlations_long['value'])
	tpr_fpr_df = pd.DataFrame({'True Positive Rate': tpr, 'False Positive Rate': fpr, 'threshold': threshold})
	roc_auc = auc(fpr, tpr)
	tpr_fpr_df["Day"]=day+", AUC="+str(round(roc_auc,2))
	return tpr_fpr_df

tpr_fpr_df_CRISPRko_d4=get_tpr_fpr_df(day4_SCEPTRE_results_CRISPRko,day="4")
tpr_fpr_df_CRISPRi_d4=get_tpr_fpr_df(day4_SCEPTRE_results_CRISPRi,day="4")
tpr_fpr_df_CRISPRko_d7=get_tpr_fpr_df(day7_SCEPTRE_results_CRISPRko,day="7")
tpr_fpr_df_CRISPRi_d7=get_tpr_fpr_df(day7_SCEPTRE_results_CRISPRi,day="7")
tpr_fpr_df_CRISPRko_d10=get_tpr_fpr_df(day10_SCEPTRE_results_CRISPRko,day="10")
tpr_fpr_df_CRISPRi_d10=get_tpr_fpr_df(day10_SCEPTRE_results_CRISPRi,day="10")
tpr_fpr_df_CRISPRko_d14=get_tpr_fpr_df(day14_SCEPTRE_results_CRISPRko,day="14")
tpr_fpr_df_CRISPRi_d14=get_tpr_fpr_df(day14_SCEPTRE_results_CRISPRi,day="14")


tpr_fpr_df_ko=pd.concat([tpr_fpr_df_CRISPRko_d14,tpr_fpr_df_CRISPRko_d10,tpr_fpr_df_CRISPRko_d7,tpr_fpr_df_CRISPRko_d4])

plt.subplots(figsize=(4, 4))
sns.lineplot(data=tpr_fpr_df_ko, x='False Positive Rate',y='True Positive Rate',
	style="Day",dashes=True,size="Day",
    errorbar=None,color="dodgerblue")
sns.despine()
plt.title("CRISPRko")
plt.savefig("../figures/CRISPRko_ROC_signedpvalcorrelation.pdf",bbox_inches="tight",dpi=600)

tpr_fpr_df_i=pd.concat([tpr_fpr_df_CRISPRi_d14,tpr_fpr_df_CRISPRi_d10,tpr_fpr_df_CRISPRi_d7,tpr_fpr_df_CRISPRi_d4])

plt.subplots(figsize=(4, 4))
sns.lineplot(data=tpr_fpr_df_i, x='False Positive Rate',y='True Positive Rate',
	style="Day",dashes=True,size="Day",
    errorbar=None,color="limegreen")
sns.despine()
plt.title("CRISPRi")
plt.savefig("../figures/CRISPRi_ROC_signedpvalcorrelation.pdf",bbox_inches="tight",dpi=600)

#repeat for A549

day7_A549_SCEPTRE_results_CRISPRko,day7_A549_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d7_A549_SCEPTRE_results.csv")

tpr_fpr_df_CRISPRko_d7_A549=get_tpr_fpr_df(day7_A549_SCEPTRE_results_CRISPRko,day="CRISPRko")
tpr_fpr_df_CRISPRi_d7_A549=get_tpr_fpr_df(day7_A549_SCEPTRE_results_CRISPRi,day="CRISPRi")

tpr_fpr_df_A549=pd.concat([tpr_fpr_df_CRISPRko_d7_A549,tpr_fpr_df_CRISPRi_d7_A549])
tpr_fpr_df_A549["Modality"]=tpr_fpr_df_A549["Day"]


plt.subplots(figsize=(4, 4))
sns.lineplot(data=tpr_fpr_df_A549, x='False Positive Rate',
             y='True Positive Rate',hue="Modality",errorbar=None,palette=["dodgerblue","limegreen"])
sns.despine()
plt.title("ROC: signed p-value correlation of guides\ntargeting the same vs. different gene\nDay 7 A549")
plt.savefig("../figures/A549_ROC_signedpvalcorrelation.pdf",bbox_inches="tight",dpi=600)


