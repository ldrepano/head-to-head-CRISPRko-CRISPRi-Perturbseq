import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

color_dictionary={"CRISPRko": 'dodgerblue', "CRISPRi": 'limegreen'}

library=pd.read_csv("../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
CRISPRi_only_guides=library[library["modality"]=="CRISPRi"]["guide_id_long"].tolist()
CRISPRko_only_guides=library[library["modality"]=="CRISPRko"]["guide_id_long"].tolist()
library["Target"]=library["guide_id_long"].apply(lambda x: x.split("_")[0])
library=library[["Target","guide_id_long","modality"]]

#absent guides receive count "0""
def get_full_library_guide_counts(sample_name):
	df=pd.read_csv("CellCountsperGuide/"+sample_name)
	full_library_counts=df.merge(library,left_on="Guide",right_on="guide_id_long",how="right")
	full_library_counts=full_library_counts[full_library_counts["modality"].isin(["CRISPRko","CRISPRi"])]
	full_library_counts=full_library_counts.fillna(0)
	return(full_library_counts)

K562_Day4=get_full_library_guide_counts("K562_Day4")
K562_Day7=get_full_library_guide_counts("K562_Day7")
K562_Day10=get_full_library_guide_counts("K562_Day10")
K562_Day14_CRISPRko=get_full_library_guide_counts("K562_Day14_CRISPRko")
K562_Day14_CRISPRko=K562_Day14_CRISPRko[K562_Day14_CRISPRko["guide_id_long"].isin(CRISPRi_only_guides)==False]
K562_Day14_CRISPRi=get_full_library_guide_counts("K562_Day14_CRISPRi")
K562_Day14_CRISPRi=K562_Day14_CRISPRi[K562_Day14_CRISPRi["guide_id_long"].isin(CRISPRko_only_guides)==False]
K562_Day14=pd.concat([K562_Day14_CRISPRko,K562_Day14_CRISPRi])
A549_Day7=get_full_library_guide_counts("A549_Day7")


def plot_counts(sample_name,df):
	plt.figure(figsize=(7,3))
	df['Target_mean'] = df.groupby('Target')['n'].transform('mean')
	df=df.sort_values(by="Target_mean",ascending=False)
	sns.stripplot(x=df["Target"],y=df["n"],hue=df["modality"],jitter=False,palette=color_dictionary)
	plt.xticks(rotation=90)
	plt.title(sample_name)
	plt.ylabel("# Cells assigned per guide")
	plt.xlabel("Guide target")
	plt.savefig("figures/"+sample_name+".png",dpi=600,bbox_inches="tight")


plot_counts("K562 Day4",K562_Day4)
plot_counts("K562 Day7",K562_Day7)
plot_counts("K562 Day10",K562_Day10)
plot_counts("K562 Day14",K562_Day14)
plot_counts("A549 Day7",A549_Day7)

def get_percent_nonesstargets_recovered(df):
	df=df[df["Target"].isin(["RPL9","DCAF13"])==False]
	return round(100*len(df[df["n"]>100])/len(df),1)

samples=["Day 4","Day 7","Day 10","Day14"]
percent_nonesstargets_recovered=(
	[get_percent_nonesstargets_recovered(K562_Day4),get_percent_nonesstargets_recovered(K562_Day7),
	get_percent_nonesstargets_recovered(K562_Day10),get_percent_nonesstargets_recovered(K562_Day14)])

percent_nonesstargets_recovered_df=pd.DataFrame({"Sample":samples,"% of library recovered\n(excluding essential-targeting)":percent_nonesstargets_recovered})

plt.figure(figsize=(5,4))
ax = sns.barplot(x=percent_nonesstargets_recovered_df["Sample"],
	y=percent_nonesstargets_recovered_df["% of library recovered\n(excluding essential-targeting)"],
	color="black")
plt.ylim(0,108) 
# Add labels to each bar
for container in ax.containers:
    ax.bar_label(container, padding=-1,fontsize=15)
plt.title("K562 Perturb-seq recovery rates")
plt.savefig("figures/K562_percentlibraryrecovered_persample.png",dpi=600,bbox_inches="tight")