#11/17/26
#plotting histograms of crossguide correlations on Day 14

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from sklearn.metrics import roc_curve, auc
import warnings
warnings.simplefilter(action='ignore', category=Warning)
plt.rc('pdf', fonttype=42)


#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]
guide_to_shortid=guideref.set_index("guide_id_long").to_dict()["guide_id_short"]


day14_SCEPTRE_results_CRISPRko=pd.read_csv("../SCEPTRE-results/d14_K562_CRISPRko_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRko["log(foldchange)"]=day14_SCEPTRE_results_CRISPRko["fold_change"].apply(lambda x: np.log10(x+0.1))
day14_SCEPTRE_results_CRISPRko["guide_id"]=day14_SCEPTRE_results_CRISPRko['grna_id'].apply(lambda x: guide_to_shortid[x])
day14_SCEPTRE_results_CRISPRko['modality']='CRISPRko'
day14_SCEPTRE_results_CRISPRko=day14_SCEPTRE_results_CRISPRko.dropna()

day14_SCEPTRE_results_CRISPRi=pd.read_csv("../SCEPTRE-results/d14_K562_CRISPRi_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRi["log(foldchange)"]=day14_SCEPTRE_results_CRISPRi["fold_change"].apply(lambda x: np.log10(x+0.1))
day14_SCEPTRE_results_CRISPRi["guide_id"]=day14_SCEPTRE_results_CRISPRi['grna_id'].apply(lambda x: guide_to_shortid[x])
day14_SCEPTRE_results_CRISPRi['modality']='CRISPRi'
day14_SCEPTRE_results_CRISPRi=day14_SCEPTRE_results_CRISPRi.dropna()

day14_SCEPTRE_results=pd.concat([day14_SCEPTRE_results_CRISPRko,day14_SCEPTRE_results_CRISPRi])

def get_signed_negative_log_pval(SCEPTRE_results_df):
	#creates a column in the supplied df called signed_logpval
	#represents each guide/transcript pair such that magnitude indicates significance of perturbation effect relative to negative controls
	#... and sign (negative vs positive) indicates upregulation vs downregulation
	SCEPTRE_results_df["effect_sign"]=SCEPTRE_results_df["fold_change"].apply(lambda x: -1 if x<1 else 1)
	SCEPTRE_results_df["-logpval"]=-np.log10(SCEPTRE_results_df["p_value"])
	SCEPTRE_results_df["signed_logpval"]=SCEPTRE_results_df["effect_sign"]*SCEPTRE_results_df["-logpval"]
	SCEPTRE_results_df=SCEPTRE_results_df.drop(["effect_sign","-logpval"],axis=1)
	return SCEPTRE_results_df

def get_correlations(SCEPTRE_results,day):
    SCEPTRE_results=get_signed_negative_log_pval(SCEPTRE_results)
    signedlogpvals=SCEPTRE_results[["response_id","grna_id","signed_logpval"]]
    signedlogpvals=pd.pivot_table(signedlogpvals, values="signed_logpval", index="response_id", columns=["grna_id"]).reset_index()
    #get pearson correlation in signed pvals attached to response gene for each guide 
    correlations = signedlogpvals.corr(numeric_only=True)
    
    correlations_long=pd.melt(correlations,ignore_index=False)
    correlations_long["other_grna_id"]=correlations_long.index
    correlations_long=correlations_long[correlations_long["grna_id"]!=correlations_long["other_grna_id"]]
    correlations_long["same_target_gene"]=correlations_long.apply(lambda x: x["grna_id"].split("_")[0]==x["other_grna_id"].split("_")[0],axis=1)

    correlations_long['Day']=day

    #removing duplicate pairs
    correlations_long['grna_id_min']=correlations_long[['grna_id','other_grna_id']].min(axis=1)
    correlations_long['grna_id_max']=correlations_long[['grna_id','other_grna_id']].max(axis=1)
    correlations_long['grna_id']=correlations_long['grna_id_min']
    correlations_long['other_grna_id']=correlations_long['grna_id_max']
    
    correlations_long=correlations_long[["grna_id","other_grna_id","same_target_gene","value"]].drop_duplicates()
    
    return correlations_long

def get_crossmodality_correlations(SCEPTRE_results):
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


	correlations_long["same_target_gene"]=correlations_long.apply(lambda x: x["grna_id"].split("_")[0]==x["guide2"].split("_")[0],axis=1)
	return correlations_long

correlations_CRISPRko_d14=get_correlations(day14_SCEPTRE_results_CRISPRko,day=14)
correlations_CRISPRi_d14=get_correlations(day14_SCEPTRE_results_CRISPRi,day=14)
correlations_crossmodality_d14=get_crossmodality_correlations(day14_SCEPTRE_results)

#plot CRISPRko
plot_df=correlations_CRISPRko_d14.copy()
plot_df['Same target gene?']=plot_df["same_target_gene"]
plt.figure()
ax=sns.histplot(data=plot_df,x='value',hue='Same target gene?',bins=100,palette={False:'lightgray',True:'blue'})
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.xlabel('Correlation between guide pairs',fontsize=12)
plt.ylabel('# guide pairs',fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Day 14 CRISPRko',fontsize=12)
plt.savefig("../figures/crisprko_correlation_hist.png",bbox_inches="tight",dpi=600)
plt.savefig("../figures/crisprko_correlation_hist.pdf",bbox_inches="tight",dpi=600)

#plot CRISPRi
plot_df=correlations_CRISPRi_d14.copy()
plot_df['Same target gene?']=plot_df["same_target_gene"]
plt.figure()
ax=sns.histplot(data=plot_df,x='value',hue='Same target gene?',bins=100,palette={False:'lightgray',True:'blue'})
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.xlabel('Correlation between guide pairs',fontsize=12)
plt.ylabel('# guide pairs',fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Day 14 CRISPRi',fontsize=12)
plt.savefig("../figures/crispri_correlation_hist.png",bbox_inches="tight",dpi=600)
plt.savefig("../figures/crispri_correlation_hist.pdf",bbox_inches="tight",dpi=600)

#plot crossmodality
plot_df=correlations_crossmodality_d14.copy()
plot_df['Same target gene?']=plot_df["same_target_gene"]
plt.figure()
ax=sns.histplot(data=plot_df,x='value',hue='Same target gene?',bins=100,palette={False:'lightgray',True:'blue'})
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.xlabel('Correlation between guide pairs',fontsize=12)
plt.ylabel('# guide pairs',fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Day 14 cross-modality',fontsize=12)
plt.savefig("../figures/crossmodality_correlation_hist.png",bbox_inches="tight",dpi=600)
plt.savefig("../figures/crossmodality_correlation_hist.pdf",bbox_inches="tight",dpi=600)