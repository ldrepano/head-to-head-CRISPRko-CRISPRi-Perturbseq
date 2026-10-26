#5/29/2026
#identifying promiscuous guides; ~=those poorly correlated with others targeting the same gene but highly reproducible

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LinearRegression
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

def processdata(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results=SCEPTRE_results.dropna()

	#creates a column in the supplied df called signed_logpval
	#represents each guide/transcript pair such that magnitude indicates significance of perturbation effect relative to negative controls
	#... and sign (negative vs positive) indicates upregulation vs downregulation
	SCEPTRE_results["effect_sign"]=SCEPTRE_results["fold_change"].apply(lambda x: -1 if x<1 else 1)
	SCEPTRE_results["-logpval"]=-np.log10(SCEPTRE_results["p_value"])
	SCEPTRE_results["signed_logpval"]=SCEPTRE_results["effect_sign"]*SCEPTRE_results["-logpval"]
	SCEPTRE_results=SCEPTRE_results.drop(["effect_sign","-logpval"],axis=1)

	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]

	#transform to one row per response gene, one column per guide
	SCEPTRE_results_CRISPRko=pd.pivot_table(SCEPTRE_results_CRISPRko[["response_id","grna_id","signed_logpval"]], values="signed_logpval", index="response_id", columns=["grna_id"])
	SCEPTRE_results_CRISPRi=pd.pivot_table(SCEPTRE_results_CRISPRi[["response_id","grna_id","signed_logpval"]], values="signed_logpval", index="response_id", columns=["grna_id"])


	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi)

day7_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d7_K562_SCEPTRE_results.csv")
day10_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d10_K562_SCEPTRE_results.csv")

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

	#repeat for timepoint2
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



d7_v_10_CRISPRko_correlations=get_correlations(day7_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRko)
d7_v_10_CRISPRi_correlations=get_correlations(day7_SCEPTRE_results_CRISPRi,day10_SCEPTRE_results_CRISPRi)


# plot modalities separately
def plot_modalities_separately(CRISPRko_df, CRISPRi_df):
	(fig,ax)=plt.subplots(1,2,figsize=(10,5),layout="constrained")
	plt.suptitle("Guide agreement with others targeting the same gene vs. reproducibility")

	modality="CRISPRko"
	ax[0].axline([0, 0], slope=1, color='black', linestyle='--',label="x=y")
	ax[0].legend(loc="lower right")
	ax[0].scatter(data=CRISPRko_df,x="median_corr",y="guide_crosstime_corr",c=color_dictionary[modality])

	ax[0].set_xlabel("Correlation with other guides targeting the same gene")
	ax[0].set_ylabel("Guide correlation Day 7 vs Day 10")
	ax[0].set_xlim(-0.01,1.01)
	ax[0].set_ylim(-0.01,1.01)
	ax[0].set_title(modality)

	modality="CRISPRi"
	ax[1].axline([0, 0], slope=1, color='black', linestyle='--',label="x=y")
	ax[1].legend(loc="lower right")
	ax[1].scatter(data=CRISPRi_df,x="median_corr",y="guide_crosstime_corr",c=color_dictionary[modality])

	ax[1].set_xlabel("Correlation with other guides targeting the same gene")
	ax[1].set_ylabel("Guide correlation Day 7 vs Day 10")
	ax[1].set_xlim(-0.01,1.01)
	ax[1].set_ylim(-0.01,1.01)
	ax[1].set_title(modality)

	plt.savefig("../figures/d7_vs_d10CRISPRkoCRISPRi_vs_withingene_correlation_K562.pdf",bbox_inches="tight",dpi=600)

plot_modalities_separately(d7_v_10_CRISPRko_correlations,d7_v_10_CRISPRi_correlations)

#plot modalities together to identify promiscuous guides (high positive residuals)
def get_scatterplot_with_outliers(df):
	df["modality_color"]=[color_dictionary[modality] for modality in df["modality"].tolist()]

	plt.figure(figsize=(5,5))
	plt.scatter(
		x=df["median_corr"],
		y=df["guide_crosstime_corr"],
		c=df["modality_color"])


	#identify outliers in plot of cross-guide vs cross sample correlation (potential promiscuous guides)
	# identify outliers with both modalities at once

	#get vertical residuals from x=y line 
	df['residual'] = (df['guide_crosstime_corr'] - df['median_corr'])

	# zscore residuals
	residual_mean=df['residual'].mean()
	residual_sd=df['residual'].std()
	df["z_scored_residuals"]=(df['residual']-residual_mean)/residual_sd

	#residual that yields a z-score of 2
	z_2=2*residual_sd+residual_mean

	#plot fitted regression and outlier cutoff 
	plt.axline(xy1=(0,0), slope=1, color='black', linestyle='-',label="x=y")
	plt.axline(xy1=(0,0+z_2), slope=1, color='black', linestyle='--',label="Residual Z-score = 2")
	plt.legend()

	# Filter for positive outliers
	positive_outliers = df[(df['z_scored_residuals'] > 2) & (df['residual'] > 0)]

	promiscuous_guide_labels=[guide_to_shortid[guide] for guide in positive_outliers["other_grna_id"].tolist()]
	promiscuous_guide_targetcorr=positive_outliers["median_corr"].tolist()
	promiscuous_guide_modalitycolor=positive_outliers["modality_color"].tolist()
	texts=[]
	promiscuous_guide_repcorr=positive_outliers["guide_crosstime_corr"].tolist()
	for i, txt in enumerate(promiscuous_guide_labels):
		texts.append(plt.text(promiscuous_guide_targetcorr[i]+0.005, promiscuous_guide_repcorr[i]+0.005, txt,c=promiscuous_guide_modalitycolor[i]))
		
	plt.xlabel("Correlation with other guides targeting the same gene")
	plt.ylabel("Correlation Day 7 vs. Day 10")
	plt.xlim(-0.01,1.01)
	plt.ylim(-0.01,1.01)
	plt.title("CRISPRko + CRISPRi")
	plt.savefig("../figures/d7_vs_d10_identifyoutliers.pdf",bbox_inches="tight",dpi=600)

	return df


d7_v_10_CRISPRko_correlations["modality"]= "CRISPRko"
d7_v_10_CRISPRi_correlations["modality"]= "CRISPRi"
d7_v_10_correlations=pd.concat([d7_v_10_CRISPRko_correlations,d7_v_10_CRISPRi_correlations])
d7_v_10_correlations_withresiduals=get_scatterplot_with_outliers(d7_v_10_correlations)

#examine association between putative off-target activity and abundance of GGs in seed 
d7_v_10_correlations_withresiduals["seed"]=d7_v_10_correlations_withresiduals["other_grna_id"].apply(lambda x: x[-12:])
d7_v_10_correlations_withresiduals["num_GG_in_seed"]=d7_v_10_correlations_withresiduals["seed"].apply(lambda x: x.count("GG"))

CRISPRidata_withresiduals=d7_v_10_correlations_withresiduals[d7_v_10_correlations_withresiduals["modality"]=="CRISPRi"]

plt.figure(figsize=(5,5))
plt.scatter(data=CRISPRidata_withresiduals[CRISPRidata_withresiduals["num_GG_in_seed"]==0],
	x="median_corr",y="guide_crosstime_corr",
	c="white",edgecolors='lightgrey', linewidths=0.5,
	s=60,label=0)
plt.scatter(data=CRISPRidata_withresiduals[CRISPRidata_withresiduals["num_GG_in_seed"]==1],
	x="median_corr",y="guide_crosstime_corr",
	c="limegreen",alpha=0.2,
	s=60,label=1)
plt.scatter(data=CRISPRidata_withresiduals[CRISPRidata_withresiduals["num_GG_in_seed"]==2],
	x="median_corr",y="guide_crosstime_corr",
	c="limegreen",s=60,label=2)
plt.scatter(data=CRISPRidata_withresiduals[CRISPRidata_withresiduals["num_GG_in_seed"]==3],
	x="median_corr",y="guide_crosstime_corr",
	c="darkgreen",s=60,label=3)
plt.legend(title="#GG in 12bp seed")
plt.axline([0, 0], slope=1, color='black', linestyle='--')
plt.xlabel("Correlation with other guides targeting the same gene")
plt.ylabel("Guide correlation Day 7 vs Day 10")
plt.xlim(-0.01,1.01)
plt.ylim(-0.01,1.01)
plt.title("Sequence composition of CRISPRi guides")
plt.savefig("../figures/CRISPRi_offtargetactivity_vsGGabundance.pdf",bbox_inches="tight",dpi=600)


def get_CRISPRi_target_kd_per_guide(sample):
    #get SCEPTRE results
    SCEPTRE_results=pd.read_csv("../SCEPTRE-results/"+sample+"_SCEPTRE_results.csv")
    SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
    SCEPTRE_results=SCEPTRE_results.dropna()
    SCEPTRE_results['sample'] = sample
    
    #add modality label, subset to CRISPRi
    SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
    SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]
    
    #get target lfc and pvalue
    SCEPTRE_results_CRISPRi["log(foldchange)"]=SCEPTRE_results_CRISPRi["fold_change"].apply(lambda x: np.log10(x+0.1))
    SCEPTRE_target_results_CRISPRi = SCEPTRE_results_CRISPRi[SCEPTRE_results_CRISPRi['response_id']==SCEPTRE_results_CRISPRi['grna_target']]
    SCEPTRE_target_results_CRISPRi=SCEPTRE_target_results_CRISPRi[['grna_id','grna_target','log(foldchange)','sample']].drop_duplicates()
    return SCEPTRE_target_results_CRISPRi

samplenames = ["d7_K562","d10_K562"]
CRISPRi_target_kd_df = pd.concat([get_CRISPRi_target_kd_per_guide(samplename) for samplename in samplenames]).reset_index(drop=True)
CRISPRi_target_kd_df["guide_id"]=CRISPRi_target_kd_df["grna_id"].apply(lambda x: guide_to_shortid[x])

CRISPRi_target_kd_with_correlations = CRISPRi_target_kd_df.merge(CRISPRidata_withresiduals,left_on='grna_id',right_on='other_grna_id')
promiscuous_CRISPRi_labels=["STK24_3","SH3KBP1_4","PYGB_2","IFT25_3","AKT2_2", "TMSB10_1","PDHA1_3"]
promiscuous_guide_targets=[guide.split("_")[0] for guide in promiscuous_CRISPRi_labels]
promiscuous_guide_targets=sorted(promiscuous_guide_targets)

plot_df = CRISPRi_target_kd_with_correlations.groupby(['guide_id','grna_target','residual'])['log(foldchange)'].mean().reset_index()

(fig,ax)=plt.subplots(2,
                      4,
                      figsize=(12,6),
                      layout="constrained")

for n,gene in enumerate(promiscuous_guide_targets):
    i=min(int(n/3),1)
    if i==0:
        j=n%3+1
    else:
        j=n-3
    df = plot_df[(plot_df['grna_target']==gene)]
    sns.scatterplot(ax=ax[i,j],
                    data=df[df['guide_id'].isin(promiscuous_CRISPRi_labels)==False],
                    x='log(foldchange)',y='residual',color='black',s=100)
    sns.scatterplot(ax=ax[i,j],
                    data=df[df['guide_id'].isin(promiscuous_CRISPRi_labels)],
                    x='log(foldchange)',y='residual',color='blue',marker="*",s=400)
    ax[i,j].set_title(gene,fontsize=16)
    ax[i,j].set_xlabel('')
    if j==0 or (i==0 and j==1):
        ax[i,j].set_ylabel('Off-target score',fontsize=16)
    else:
        ax[i,j].set_ylabel('')
    ax[i,j].tick_params(axis="both", labelsize=12)

ax[0, 0].set_visible(False)
fig.supxlabel('Target LFC (average across Day 7 and Day 10)',fontsize=16)

plt.savefig("../figures/target_lfc_per_guide_per_target_with_promiscuous_guide_vs_offtarget_score.pdf",bbox_inches="tight",dpi=600)