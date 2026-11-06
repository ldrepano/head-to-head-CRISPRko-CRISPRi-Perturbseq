#10/2/2026 (revised: rank by signed p-value, fold-change tie-break, ECDF)
#Day 4,7,10,14 K562 examine knockdown of target gene relative to intergenic controls across CRISPRko and CRISPRi,
#percentile of the target's (signed) p-value compared to all other transcripts tested against the same guide

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=Warning)
plt.rc('pdf', fonttype=42)

P_FLOOR = 1e-300          #avoid log10(0) if any p-values are reported as exactly 0

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
	SCEPTRE_results=SCEPTRE_results.dropna(subset=["p_value","fold_change"])
	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"].copy()
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"].copy()
	return (SCEPTRE_results_CRISPRko,SCEPTRE_results_CRISPRi)

day4_SCEPTRE_results_CRISPRko,day4_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d4_K562_SCEPTRE_results.csv")
day7_SCEPTRE_results_CRISPRko,day7_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d7_K562_SCEPTRE_results.csv")
day10_SCEPTRE_results_CRISPRko,day10_SCEPTRE_results_CRISPRi=processdata_nonday14("../SCEPTRE-results/d10_K562_SCEPTRE_results.csv")

#load/process day 14 data
def processdata(sceptre_results_filepath):
	SCEPTRE_results=pd.read_csv(sceptre_results_filepath)
	SCEPTRE_results["log(foldchange)"]=SCEPTRE_results["fold_change"].apply(lambda x: np.log10(x+0.1))
	SCEPTRE_results=SCEPTRE_results.dropna(subset=["p_value","fold_change"])
	return SCEPTRE_results

day14_SCEPTRE_results_CRISPRko=processdata("../SCEPTRE-results/d14_K562_CRISPRko_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRi=processdata("../SCEPTRE-results/d14_K562_CRISPRi_SCEPTRE_results.csv")
day14_SCEPTRE_results_CRISPRko["modality"]="CRISPRko"
day14_SCEPTRE_results_CRISPRi["modality"]="CRISPRi"

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

# lower ranking = stronger evidence of knockdown 
logp = np.log10(modalities_combined['p_value'].clip(lower=P_FLOOR))
modalities_combined['rank_score'] = np.where(modalities_combined['fold_change'] < 1, logp, -logp)

#percentile rank of the target within each guide: primary key = rank_score, tie-break = fold change
#combined via dense ranks so that identical (score, fold change) pairs still tie with method='min'
grp = ['modality', 'grna_id', 'Day']
r_score = modalities_combined.groupby(grp)['rank_score'].rank(method='dense')
r_fc = modalities_combined.groupby(grp)['fold_change'].rank(method='dense')
modalities_combined['rank_key'] = r_score * (r_fc.max() + 1) + r_fc # weights p-val most and then ties can be broken by fc 
modalities_combined['ptl'] = modalities_combined.groupby(grp)['rank_key'].rank(method='min', pct=True) # min rank method 

#filter down to only the target transcripts
target_modalities_combined = modalities_combined[modalities_combined["response_id"] == modalities_combined["grna_target"]].copy()
target_modalities_combined['ptl_100'] = target_modalities_combined['ptl'] * 100

#tie report: how many transcripts (including the target) share the target's ranking key for that guide
#n_tied_p   = ties on p-value alone (before the fold-change tie-break)
#n_tied     = ties remaining after the tie-break (these are what actually affect the percentile)
tie_counts_p = (modalities_combined
    .groupby(grp + ['rank_score']).size().rename('n_tied_p').reset_index())
tie_counts = (modalities_combined
    .groupby(grp + ['rank_key']).size().rename('n_tied').reset_index())
target_modalities_combined = (target_modalities_combined
    .merge(tie_counts_p, on=grp + ['rank_score'], how='left')
    .merge(tie_counts, on=grp + ['rank_key'], how='left'))

tie_summary = (target_modalities_combined
    .groupby(['modality', 'Day'])
    .agg(n_guides=('grna_id', 'size'),
         frac_target_upregulated=('fold_change', lambda s: (s >= 1).mean()),
         frac_target_tied_on_p=('n_tied_p', lambda s: (s > 1).mean()),
         frac_target_tied_after_fc=('n_tied', lambda s: (s > 1).mean()),
         median_tie_size_when_tied=('n_tied', lambda s: s[s > 1].median() if (s > 1).any() else np.nan)))
print("Tie report (target transcript):")
print(tie_summary)

target_modalities_combined = target_modalities_combined.sort_values(by='ptl_100')
target_modalities_combined_plot = target_modalities_combined[['ptl_100','modality','Day']]

#empirical CDF: y = fraction of guides with ptl_100 <= x (ties handled correctly)
#curve is anchored at (0, 0) and extended to (100, 1) for visualization
def ecdf_xy(values, xmin=0, xmax=100):
    x = np.sort(np.asarray(values))
    y = np.arange(1, len(x) + 1) / len(x)
    x = np.concatenate([[xmin], x, [xmax]])
    y = np.concatenate([[0], y, [1]])
    return x, y

day_order = ["14", "10", "7", "4"]

def plot_cum_dist(df, color, title, outpath):
    #build ECDF points per day
    curves = []
    for day in day_order:
        vals = df.loc[df['Day'] == day, 'ptl_100']
        if len(vals) == 0:
            continue
        x, y = ecdf_xy(vals)
        curves.append(pd.DataFrame({'ptl_100': x, 'frac': y, 'Day': day}))
    curves = pd.concat(curves, ignore_index=True)
    days_present = [d for d in day_order if d in set(curves['Day'])]

    fig, ax = plt.subplots(figsize=(4, 4))
    sns.lineplot(data=curves, x='ptl_100', y='frac',
                 style='Day', size='Day', dashes=True,
                 style_order=days_present, size_order=days_present,
                 estimator=None, sort=False,      #no averaging of ties, keep ECDF order
                 drawstyle='steps-post',          #step function
                 color=color, ax=ax)
    ax.axvline(x=20, color='red', linestyle='--')
    ax.set_xlabel("Target percentile by signed p-value\n(relative to other transcripts, same guide)")
    ax.set_ylabel("Fraction of guides at or below")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower right", title="Day")
    sns.despine()
    ax.set_title(title)
    fig.savefig(outpath, bbox_inches="tight", dpi=600)

#cumulative distribution plots
plot_cum_dist(target_modalities_combined_plot.loc[target_modalities_combined_plot['modality']=="CRISPRko"],
              "dodgerblue", "CRISPRko", "../figures/pval_CRISPRko_cum_dist_target_percentile_ecdf.pdf")
plot_cum_dist(target_modalities_combined_plot.loc[target_modalities_combined_plot['modality']=="CRISPRi"],
              "limegreen", "CRISPRi", "../figures/pval_CRISPRi_cum_dist_target_percentile_ecdf.pdf")

#fraction of guides with target in bottom 20% (the red line), for reference
print("\nFraction of guides with target percentile <= 20:")
print(target_modalities_combined.groupby(['modality', 'Day'])['ptl_100'].apply(lambda s: (s <= 20).mean()))

#sanity check: best-ranked CRISPRi Day 7 targets
d7i = target_modalities_combined.query("modality=='CRISPRi' and Day=='7'").sort_values('ptl_100')
print(d7i[['grna_id','response_id','p_value','fold_change','ptl_100','n_tied_p','n_tied']].head())

#how many transcripts share the minimum p-value for that guide
g = d7i.iloc[0]['grna_id']
sub = modalities_combined.query("modality=='CRISPRi' and Day=='7' and grna_id==@g")
print((sub['p_value'] == sub['p_value'].min()).sum(), len(sub))