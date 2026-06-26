#5/29/26
#identifying log fold change of target knockdown from SCEPTRE results relative to all transcripts

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
plt.rc('pdf', fonttype=42)

import warnings
warnings.simplefilter(action='ignore', category=Warning)


#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]

#cannot distinguish these since CRISPRko, CRISPRi perturbed cells are sequenced together
guides_repeated_across_libraries=["COTL1_TACAACCTGGTGCGCGACGA","MOSMO_GGCGATGGCGAAGATATCGG"]

def get_frac_kd(sample):
	SCEPTRE_results=pd.read_csv("../SCEPTRE-results/d"+sample+"_K562_SCEPTRE_results.csv")
	SCEPTRE_results=SCEPTRE_results[SCEPTRE_results["grna_id"].isin(guides_repeated_across_libraries)==False]
	SCEPTRE_results=SCEPTRE_results.dropna()
	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	SCEPTRE_results_CRISPRko=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRko"]
	SCEPTRE_results_CRISPRi=SCEPTRE_results[SCEPTRE_results["modality"]=="CRISPRi"]


	SCEPTRE_results_CRISPRko_target=SCEPTRE_results_CRISPRko[SCEPTRE_results_CRISPRko["response_id"]==SCEPTRE_results_CRISPRko["grna_target"]].reset_index()
	CRISPRko_frac_target_kd= len(SCEPTRE_results_CRISPRko_target[SCEPTRE_results_CRISPRko_target["significant"]])/len(SCEPTRE_results_CRISPRko_target)

	SCEPTRE_results_CRISPRi_target=SCEPTRE_results_CRISPRi[SCEPTRE_results_CRISPRi["response_id"]==SCEPTRE_results_CRISPRi["grna_target"]].reset_index()
	CRISPRi_frac_target_kd= len(SCEPTRE_results_CRISPRi_target[SCEPTRE_results_CRISPRi_target["significant"]])/len(SCEPTRE_results_CRISPRi_target)

	return (CRISPRko_frac_target_kd,CRISPRi_frac_target_kd)


CRISPRko_frac_target_kd_list=[]
CRISPRi_frac_target_kd_list=[]
samples=["4","7","10"]
for sample in samples:
	result=get_frac_kd(sample)
	CRISPRko_frac_target_kd_list.append(100*result[0])
	CRISPRi_frac_target_kd_list.append(100*result[1])

#add day14 data
samples.append("14")

SCEPTRE_results_d14CRISPRko=pd.read_csv("../SCEPTRE-results/d14_K562_CRISPRko_SCEPTRE_results.csv")
SCEPTRE_results_d14CRISPRko=SCEPTRE_results_d14CRISPRko.dropna()
SCEPTRE_results_d14CRISPRko_target=SCEPTRE_results_d14CRISPRko[SCEPTRE_results_d14CRISPRko["response_id"]==SCEPTRE_results_d14CRISPRko["grna_target"]].reset_index()
d14CRISPRko_frac_target_kd= len(SCEPTRE_results_d14CRISPRko_target[SCEPTRE_results_d14CRISPRko_target["significant"]])/len(SCEPTRE_results_d14CRISPRko_target)
CRISPRko_frac_target_kd_list.append(100*d14CRISPRko_frac_target_kd)

SCEPTRE_results_d14CRISPRi=pd.read_csv("../SCEPTRE-results/d14_K562_CRISPRi_SCEPTRE_results.csv")
SCEPTRE_results_d14CRISPRi=SCEPTRE_results_d14CRISPRi.dropna()
SCEPTRE_results_d14CRISPRi_target=SCEPTRE_results_d14CRISPRi[SCEPTRE_results_d14CRISPRi["response_id"]==SCEPTRE_results_d14CRISPRi["grna_target"]].reset_index()
d14CRISPRi_frac_target_kd= len(SCEPTRE_results_d14CRISPRi_target[SCEPTRE_results_d14CRISPRi_target["significant"]])/len(SCEPTRE_results_d14CRISPRi_target)
CRISPRi_frac_target_kd_list.append(100*d14CRISPRi_frac_target_kd)


frac_kd_df=pd.DataFrame({"CRISPRko":CRISPRko_frac_target_kd_list,"CRISPRi":CRISPRi_frac_target_kd_list},index=samples).reset_index()
frac_kd_df_long=frac_kd_df.melt(id_vars="index",value_name="percent_kd",var_name="Modality")
plt.figure()

ax=sns.barplot(data=frac_kd_df_long,x="index",y="percent_kd",hue="Modality",palette=["dodgerblue","limegreen"])
for i in range(len(ax.containers)):
	container=ax.containers[i]
	label=CRISPRko_frac_target_kd_list[i]
	ax.bar_label(container, fmt="%.1f", padding=0,label=label)
plt.xlabel("Day")
plt.ylabel("% guides with target knockdown")
plt.title("Target knockdown rates")
plt.ylim(0,105)
plt.legend(loc="lower left")

plt.savefig("../figures/frac_guides_with_target_kd.pdf",bbox_inches="tight",dpi=600)


