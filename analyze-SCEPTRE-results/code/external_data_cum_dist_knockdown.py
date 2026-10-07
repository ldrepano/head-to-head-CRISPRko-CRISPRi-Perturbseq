#10/2/2026
#examine knockdown of target gene relative to intergenic controls for Yao et al. CRISPRko and CRISPRi data

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import warnings
warnings.simplefilter(action='ignore', category=Warning)
plt.rc('pdf', fonttype=42)

def process_Yao_SCEPTRE_results(modality):
    SCEPTRE_results = pd.read_csv(f"../SCEPTRE-results/{modality}_Yao_SCEPTRE_results.csv")
    SCEPTRE_results["log(foldchange)"] = np.log10(SCEPTRE_results["fold_change"] + 0.1)
    SCEPTRE_results=SCEPTRE_results.dropna()
    SCEPTRE_results=SCEPTRE_results.sort_values(by="grna_id")
    
    SCEPTRE_results_target=SCEPTRE_results[SCEPTRE_results["response_id"]==SCEPTRE_results["grna_target"]].reset_index()
    guides_with_target_detected=SCEPTRE_results_target["grna_id"].unique() 
    SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_with_target_detected)].reset_index(drop=True)
    SCEPTRE_results["modality"] = SCEPTRE_results['grna_id'].apply(lambda x: x.split('_')[-1])

    
    return SCEPTRE_results,SCEPTRE_results_target

Yao_SCEPTRE_results_CRISPRko,Yao_SCEPTRE_results_CRISPRko_target=process_Yao_SCEPTRE_results("CRISPRko")
Yao_SCEPTRE_results_CRISPRi,Yao_SCEPTRE_results_CRISPRi_target=process_Yao_SCEPTRE_results("CRISPRi")

yao_combined = pd.concat([Yao_SCEPTRE_results_CRISPRko, Yao_SCEPTRE_results_CRISPRi], ignore_index=True)
yao_combined['ptl'] = yao_combined.groupby(['modality', 'grna_id'])['log(foldchange)'].rank(pct=True)

# Filter down to only the target transcripts 
yao_target_combined = yao_combined[yao_combined["response_id"] == yao_combined["grna_target"]].copy()
yao_target_combined['ptl_100'] = yao_target_combined['ptl'] * 100

yao_target_combined = yao_target_combined.sort_values(by='ptl_100')
yao_target_combined['cum_ptl_ptl'] = yao_target_combined.groupby('modality')['ptl_100'].rank(method='max', pct=True)

#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]

#cannot distinguish these at Day 4,7,10 since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

def processdata_nonday14(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log10(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna()
	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]
	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi)

day10_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d10_K562_SCEPTRE_results.csv")

# Calculate percentile rank of the fold change for each guide
day10_combined = pd.concat([day10_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRi], ignore_index = True)
day10_combined['ptl'] = day10_combined.groupby(['modality', 'grna_id'])['log(foldchange)'].rank(pct=True)

# Filter down to only the target transcripts 
target_day10_combined = day10_combined[day10_combined["response_id"] == day10_combined["grna_target"]].copy()
target_day10_combined['ptl_100'] = target_day10_combined['ptl'] * 100

target_day10_combined = target_day10_combined.sort_values(by='ptl_100')
target_day10_combined['cum_ptl_ptl'] = target_day10_combined.groupby(['modality'])['ptl_100'].rank(method='max', pct=True)

def get_frac_kd(SCEPTRE_result, modality):
    """Fraction of guides where the target transcript is called significant."""
    df = SCEPTRE_result
    target = df[df["response_id"] == df["grna_target"]].reset_index()
    if len(target) == 0:
        return np.nan
    return len(target[target["significant"]]) / len(target)

Yao_CRISPRko_pct = 100 * get_frac_kd(Yao_SCEPTRE_results_CRISPRko, "CRISPRko")
Yao_CRISPRi_pct = 100 * get_frac_kd(Yao_SCEPTRE_results_CRISPRi,"CRISPRi")
day10_CRISPRko_pct = 100 * get_frac_kd(day10_SCEPTRE_results_CRISPRko, "CRISPRko")
day10_CRISPRi_pct = 100 * get_frac_kd(day10_SCEPTRE_results_CRISPRi,"CRISPRi")

modality_colors = {
    "CRISPRko": "dodgerblue",
    "CRISPRi": "limegreen"} 
plt.subplots(figsize=(4, 4))
sns.lineplot(data=yao_target_combined, 
             x='ptl_100',
             y='cum_ptl_ptl',
             dashes=True,
             errorbar=None,
             palette =modality_colors,
             hue = 'modality')
plt.axvline(x=20, color ='red', linestyle ='--')
plt.xlabel("Target percentile\n(relative to other transcripts, same guide)")
plt.ylabel("Fraction of guides at or below")
plt.xlim(0,100)
plt.legend(loc="lower right")

sns.despine()
plt.title("Yao et al.")
plt.savefig("../figures/yao_cum_dist_target_ko_percentile.pdf",bbox_inches="tight",dpi=600)
plt.show()

datasets = ["This study, Day 10 K562", "Yao et al."]

ko_pcts = [day10_CRISPRko_pct, Yao_CRISPRko_pct]
i_pcts = [day10_CRISPRi_pct, Yao_CRISPRi_pct]

x = np.arange(len(datasets))  
width = 0.35  

fig, ax_b = plt.subplots(figsize=(4, 4))

bars1 = ax_b.bar(x - width/2, ko_pcts, width, label='CRISPRko', color=modality_colors["CRISPRko"])
bars2 = ax_b.bar(x + width/2, i_pcts, width, label='CRISPRi', color=modality_colors["CRISPRi"])

for bars in [bars1, bars2]:
    for bar in bars:
        height = float(bar.get_height())
        ax_b.text(
            bar.get_x() + bar.get_width() / 2, 
            height + 1, 
            f"{height:.1f}",
            ha="center", 
            va="bottom"
        )

ax_b.set_xticks(x)
ax_b.set_xticklabels(datasets)
ax_b.set_ylim(0, 110)  
ax_b.set_ylabel("% guides with target knockdown")
ax_b.set_title("Target knockdown rate (FDR < 0.1)")
ax_b.legend(loc="upper right")

plt.savefig("../figures/target_knockdown_rate.pdf",bbox_inches="tight",dpi=600)
plt.show()
