#6/1/2026
#identifies CRISPRko bias in hits on target chromosome arm AKT2-targeting guides

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.simplefilter(action='ignore', category=Warning)

library=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
library=library[library["guide_id_short"].str.count("INTERGENIC")==0]
library["target"]=library["guide_id_short"].apply(lambda x:x.split("_")[0])

#get all chromosome centromere positions (source info: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cytoBand)
arm_coordinates=pd.read_csv("../../reference/UCSC_hg38_chromosome_band_track.csv")
arm_coordinates=arm_coordinates.dropna()
#data gives band coordinates. get centromere coordinates
arm_coordinates["arm"]=arm_coordinates["name"].apply(lambda x: x[0])
p_arm_coordinates=arm_coordinates[arm_coordinates["arm"]=="p"].reset_index(drop=True)
centromere_coords=p_arm_coordinates.groupby("#chrom").agg(centromere_position=("chromEnd","max")).reset_index()
centromere_coords["#chrom"]=centromere_coords["#chrom"].str.replace("chr","")
centromere_coords_map=centromere_coords.set_index("#chrom")[["centromere_position"]].to_dict()["centromere_position"]

#ensembl 115 gene locations from BioMart
gene_locs = pd.read_table('../../reference/gene_locations_ensembl115.txt', sep = '\t', names = ['gene', 'chrom', 'start', 'end'])
all_chrs=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
gene_locs=gene_locs[gene_locs["chrom"].isin(all_chrs)]
gene_locs["gene_centromere_position"]=gene_locs["chrom"].apply(lambda x: centromere_coords_map[x])
gene_locs["target_gene_arm"]=gene_locs.apply(lambda x: "p" if x["start"]<x["gene_centromere_position"] else "q",axis=1)
gene_to_chrom = gene_locs.set_index('gene')['chrom']
gene_to_pos = gene_locs.set_index('gene')['start']


#manually fixing MOSMO and MGAT4B locations -- annotated on patches for some reason
gene_to_chrom["MOSMO"]="16"
gene_to_pos["MOSMO"]= 22007638
gene_to_chrom["MGAT4B"]="5"
gene_to_pos["MGAT4B"]= 179224597

gene_to_arm=gene_locs.set_index('gene')['target_gene_arm']
gene_to_arm["MOSMO"]= "p"
gene_to_arm["MGAT4B"]="q"

library_targets=library["target"].unique().tolist()

(fix,ax)=plt.subplots(5,5,figsize=(10,5), layout="constrained")

plt.suptitle("Position along chromosome defined as 'towards telomere'\nPer library target")

for i in range(len(library_targets)):
	target=library_targets[i]
	chromosome=gene_to_chrom[target]
	row=i//5
	col=i%5

	ax[row,col].set_xlabel("Chromosome "+chromosome,style='italic')
	chr_end=gene_locs[gene_locs["chrom"]==chromosome]["end"].max()
	ax[row,col].set_xlim(0,chr_end)
	ax[row,col].set_yticks([])
	ax[row,col].set_xticks([])
	ax[row,col].set_title("")
	ax[row,col].axvline(x=gene_to_pos[target],color="blue",label=target)
	ax[row,col].legend()

	if gene_to_arm[target]=="p":
		ax[row,col].axvspan(0, gene_to_pos[target], edgecolor="grey", fill=False,hatch="////")
	else:
		ax[row,col].axvspan(gene_to_pos[target],chr_end, edgecolor="grey", fill=False,hatch="////")


plt.savefig("../figures/region_towards_telomere_per_library_target.png",bbox_inches="tight",dpi=600)







