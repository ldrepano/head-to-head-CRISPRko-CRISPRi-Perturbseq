#10/2/2026
#Day 4,7,10,14 K562 examine knockdown of target gene relative to intergenic controls across CRISPRko and CRISPRi, percentile compared to all other transcripts

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import warnings
warnings.simplefilter(action='ignore', category=Warning)
plt.rc('pdf', fonttype=42)

#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]

#cannot distinguish these at Day 4,7,10 since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

#load/process day 4,7,10 data
def processdata_nonday14(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log10(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna()
	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]
	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi)

day4_SCEPTRE_results_CRISPRko,day4_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d4_K562_SCEPTRE_results.csv")
day7_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d7_K562_SCEPTRE_results.csv")
day10_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d10_K562_SCEPTRE_results.csv")

#load/process day 14 dara
def processdata(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log10(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna()
	return SCEPTRE_results

day14_SCEPTRE_results_CRISPRko=processdata("../SCEPTRE-results/d14_K562_CRISPRko_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d14_K562_CRISPRi_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRko["modality"]="CRISPRko"
day14_SCEPTRE_results_CRISPRi["modality"]="CRISPRi"
day14_SCEPTRE_results=pd.concat([day14_SCEPTRE_results_CRISPRko,day14_SCEPTRE_results_CRISPRi])

#add column identifying day of experiment 
def add_day_col(target_dict):
    for day_label, df_list in target_dict.items():
        for df in df_list:
            df["Day"] = day_label

dayday_dict = {
    "4": [day4_SCEPTRE_results_CRISPRko, day4_SCEPTRE_results_CRISPRi],
    "7": [day7_SCEPTRE_results_CRISPRko, day7_SCEPTRE_results_CRISPRi],
    "10": [day10_SCEPTRE_results_CRISPRko, day10_SCEPTRE_results_CRISPRi],
    "14": [day14_SCEPTRE_results_CRISPRko, day14_SCEPTRE_results_CRISPRi],
}
add_day_col(dayday_dict)

#combine dataframes
CRISPRko_combined = pd.concat([day4_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRko,day14_SCEPTRE_results_CRISPRko], ignore_index=True)
CRISPRi_combined = pd.concat([day4_SCEPTRE_results_CRISPRi,day7_SCEPTRE_results_CRISPRi,day10_SCEPTRE_results_CRISPRi,day14_SCEPTRE_results_CRISPRi], ignore_index=True)
modalities_combined = pd.concat([CRISPRi_combined, CRISPRko_combined], ignore_index=True)

#calculate percentile rank of the fold change for each guide
modalities_combined['ptl'] = modalities_combined.groupby(['modality', 'grna_id', 'Day'])['log(foldchange)'].rank(pct=True)

#filter down to only the target transcripts 
target_modalities_combined = modalities_combined[modalities_combined["response_id"] == modalities_combined["grna_target"]].copy()
target_modalities_combined['ptl_100'] = target_modalities_combined['ptl'] * 100
target_modalities_combined = target_modalities_combined.sort_values(by='ptl_100')
#every guide's ptl_100 value gets ranked against all the other guides' ptl_100 values in that same group, and that rank is expressed as a percentile (pct=True)
target_modalities_combined['cum_ptl_ptl'] = target_modalities_combined.groupby(['modality', 'Day'])['ptl_100'].rank(method='max', pct=True)

target_modalities_combined_plot = target_modalities_combined[['ptl_100','modality','Day','cum_ptl_ptl']]

# add data points extending to x = 100, y = 1.0 solely for visualization purposes
new_rows_data = [
    {'ptl_100': 100, 'modality': 'CRISPRi','Day': '4','cum_ptl_ptl': 1},
    {'ptl_100': 100, 'modality': 'CRISPRi', 'Day': '7','cum_ptl_ptl': 1},
    {'ptl_100': 100, 'modality': 'CRISPRi','Day': '10','cum_ptl_ptl': 1},
    {'ptl_100': 100, 'modality': 'CRISPRi','Day': '14','cum_ptl_ptl': 1},
    {'ptl_100': 100, 'modality': 'CRISPRko', 'Day': '4','cum_ptl_ptl': 1},
    {'ptl_100': 100, 'modality': 'CRISPRko','Day': '7','cum_ptl_ptl': 1},
    {'ptl_100': 100, 'modality': 'CRISPRko', 'Day': '10','cum_ptl_ptl': 1},
    {'ptl_100': 100, 'modality': 'CRISPRko','Day': '14','cum_ptl_ptl': 1}
]
new_rows_df = pd.DataFrame(new_rows_data)
target_modalities_combined_plot = pd.concat([target_modalities_combined_plot,new_rows_df], ignore_index=True)

#cumulative distribution plots
plt.subplots(figsize=(4, 4))
sns.lineplot(data=target_modalities_combined_plot.loc[target_modalities_combined_plot['modality']=="CRISPRko"], x='ptl_100',y='cum_ptl_ptl',
	style="Day",dashes=True,size="Day",
    errorbar=None,color="dodgerblue")
plt.axvline(x=20, color ='red', linestyle ='--')
plt.xlabel("Target percentile\n(relative to other transcripts, same guide)")
plt.ylabel("Fraction of guides at or below")
plt.xlim(0,100)
plt.legend(loc="lower right", title = "Day")
sns.despine()
plt.title("CRISPRko")
plt.savefig("../figures/CRISPRko_cum_dist_target_ko_percentile.pdf",bbox_inches="tight",dpi=600)

plt.subplots(figsize=(4, 4))
sns.lineplot(data=target_modalities_combined_plot.loc[target_modalities_combined_plot['modality']=="CRISPRi"], x='ptl_100',y='cum_ptl_ptl',
	style="Day",dashes=True,size="Day",
    errorbar=None,color="limegreen")
plt.axvline(x=20, color ='red', linestyle ='--')
plt.xlabel("Target percentile\n(relative to other transcripts, same guide)")
plt.ylabel("Fraction of guides at or below")
plt.xlim(0,100)
plt.legend(loc="lower right", title = "Day")
sns.despine()
plt.title("CRISPRi")
plt.savefig("../figures/CRISPRi_cum_dist_target_ko_percentile.pdf",bbox_inches="tight",dpi=600)