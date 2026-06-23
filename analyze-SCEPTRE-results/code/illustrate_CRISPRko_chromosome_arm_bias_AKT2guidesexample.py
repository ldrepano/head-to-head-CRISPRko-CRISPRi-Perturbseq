#6/1/2026
#identifies CRISPRko bias in hits on target chromosome arm AKT2-targeting guides

import pandas as pd
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=Warning)

#first read in reference such that CRISPRko vs CRISPRi guides can be distinguished
guideref=pd.read_csv("../../reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")

guide_to_modality=guideref.set_index("guide_id_long").to_dict()["modality"]
guide_to_shortid=guideref.set_index("guide_id_long").to_dict()["guide_id_short"]

#ensembl 115 gene locations from BioMart
gene_locs = pd.read_table('../../reference/gene_locations_ensembl115.txt', sep = '\t', names = ['gene', 'chrom', 'start', 'end'])

#identify significantly downregulated genes
def get_sigdownreg_genes(filepath):
	SCEPTRE_results=pd.read_csv(filepath)
	SCEPTRE_results=SCEPTRE_results.dropna()
	SCEPTRE_results["grna_id_short"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_shortid[x])
	SCEPTRE_results=SCEPTRE_results.sort_values(by="grna_id_short")
	SCEPTRE_results["modality"]=SCEPTRE_results['grna_id'].apply(lambda x: guide_to_modality[x])
	#restrict to significant hits (FDR<0.1)
	SCEPTRE_sig=SCEPTRE_results[SCEPTRE_results["significant"]]
	#remove target from responses
	SCEPTRE_sig=SCEPTRE_sig[SCEPTRE_sig["grna_target"]!=SCEPTRE_sig["response_id"]]
	#subset to downreg
	SCEPTRE_downreg=SCEPTRE_sig[SCEPTRE_sig["fold_change"]<1]
	return SCEPTRE_downreg


day7_with_chr_info=get_sigdownreg_genes("../SCEPTRE-results/d7_K562_SCEPTRE_results.csv")
Day7_sig_CRISPRko=day7_with_chr_info[day7_with_chr_info["modality"]=="CRISPRko"]
Day7_sig_CRISPRi=day7_with_chr_info[day7_with_chr_info["modality"]=="CRISPRi"]


# -- Making Plots --


#AKT2 rugplots
AKT2_targeting_CRISPRko=Day7_sig_CRISPRko[Day7_sig_CRISPRko["grna_target"]=="AKT2"]
AKT2_targeting_CRISPRi=Day7_sig_CRISPRi[Day7_sig_CRISPRi["grna_target"]=="AKT2"]

(fix,ax)=plt.subplots(4,2,figsize=(12,4), layout="constrained")

plt.suptitle("Genes along chromosome 19 downregulated (SCEPTRE FDR<0.1) at Day 7\nAKT2 in blue")

AKT2_1_DE_genes_CRISPRko=AKT2_targeting_CRISPRko[AKT2_targeting_CRISPRko["grna_id_short"]=="AKT2_1"]["response_id"].tolist()
AKT2_1_DE_genes_CRISPRko_positions=gene_locs[gene_locs["gene"].isin(AKT2_1_DE_genes_CRISPRko)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_1_DE_genes_CRISPRko_positions:
	ax[0,0].axvline(x=gene_pos,color="grey")
ax[0,0].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[0,0].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_xlabel("")
ax[0,0].set_title("CRISPRko guide AKT2_1")
chromosome_19_length= 58616341 #according to UCSC browser hg38 
ax[0,0].set_xlim(0,chromosome_19_length)

AKT2_2_DE_genes_CRISPRko=AKT2_targeting_CRISPRko[AKT2_targeting_CRISPRko["grna_id_short"]=="AKT2_2"]["response_id"].tolist()
AKT2_2_DE_genes_CRISPRko_positions=gene_locs[gene_locs["gene"].isin(AKT2_2_DE_genes_CRISPRko)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_2_DE_genes_CRISPRko_positions:
	ax[1,0].axvline(x=gene_pos,color="grey")
ax[1,0].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[1,0].set_yticks([])
ax[1,0].set_xticks([])
ax[1,0].set_xlabel("")
ax[1,0].set_title("CRISPRko guide AKT2_2")
ax[1,0].set_xlim(0,chromosome_19_length)

AKT2_3_DE_genes_CRISPRko=AKT2_targeting_CRISPRko[AKT2_targeting_CRISPRko["grna_id_short"]=="AKT2_3"]["response_id"].tolist()
AKT2_3_DE_genes_CRISPRko_positions=gene_locs[gene_locs["gene"].isin(AKT2_3_DE_genes_CRISPRko)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_3_DE_genes_CRISPRko_positions:
	ax[2,0].axvline(x=gene_pos,color="grey")
ax[2,0].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[2,0].set_yticks([])
ax[2,0].set_xticks([])
ax[2,0].set_xlabel("")
ax[2,0].set_title("CRISPRko guide AKT2_3")
ax[2,0].set_xlim(0,chromosome_19_length)

AKT2_4_DE_genes_CRISPRko=AKT2_targeting_CRISPRko[AKT2_targeting_CRISPRko["grna_id_short"]=="AKT2_4"]["response_id"].tolist()
AKT2_4_DE_genes_CRISPRko_positions=gene_locs[gene_locs["gene"].isin(AKT2_4_DE_genes_CRISPRko)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_4_DE_genes_CRISPRko_positions:
	ax[3,0].axvline(x=gene_pos,color="grey")
ax[3,0].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[3,0].set_yticks([])
ax[3,0].set_xlabel("Position (bp) along chr19 (hg38)")
ax[3,0].set_title("CRISPRko guide AKT2_4")
ax[3,0].set_xlim(0,chromosome_19_length)


AKT2_1_DE_genes_CRISPRi=AKT2_targeting_CRISPRi[AKT2_targeting_CRISPRi["grna_id_short"]=="AKT2_1"]["response_id"].tolist()
AKT2_1_DE_genes_CRISPRi_positions=gene_locs[gene_locs["gene"].isin(AKT2_1_DE_genes_CRISPRi)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_1_DE_genes_CRISPRi_positions:
	ax[0,1].axvline(x=gene_pos,color="grey")
ax[0,1].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[0,1].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_xlabel("")
ax[0,1].set_title("CRISPRi guide AKT2_1")
ax[0,1].set_xlim(1,chromosome_19_length)

AKT2_2_DE_genes_CRISPRi=AKT2_targeting_CRISPRi[AKT2_targeting_CRISPRi["grna_id_short"]=="AKT2_2"]["response_id"].tolist()
AKT2_2_DE_genes_CRISPRi_positions=gene_locs[gene_locs["gene"].isin(AKT2_2_DE_genes_CRISPRi)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_2_DE_genes_CRISPRi_positions:
	ax[1,1].axvline(x=gene_pos,color="grey")
ax[1,1].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[1,1].set_yticks([])
ax[1,1].set_xticks([])
ax[1,1].set_xlabel("")
ax[1,1].set_title("CRISPRi guide AKT2_2")
ax[1,1].set_xlim(1,chromosome_19_length)

AKT2_3_DE_genes_CRISPRi=AKT2_targeting_CRISPRi[AKT2_targeting_CRISPRi["grna_id_short"]=="AKT2_3"]["response_id"].tolist()
AKT2_3_DE_genes_CRISPRi_positions=gene_locs[gene_locs["gene"].isin(AKT2_3_DE_genes_CRISPRi)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_3_DE_genes_CRISPRi_positions:
	ax[2,1].axvline(x=gene_pos,color="grey")
ax[2,1].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[2,1].set_yticks([])
ax[2,1].set_xticks([])
ax[2,1].set_xlabel("")
ax[2,1].set_title("CRISPRi guide AKT2_3")
ax[2,1].set_xlim(1,chromosome_19_length)

AKT2_4_DE_genes_CRISPRi=AKT2_targeting_CRISPRi[AKT2_targeting_CRISPRi["grna_id_short"]=="AKT2_4"]["response_id"].tolist()
AKT2_4_DE_genes_CRISPRi_positions=gene_locs[gene_locs["gene"].isin(AKT2_4_DE_genes_CRISPRi)&(gene_locs["chrom"]=="19")]["start"].tolist()
for gene_pos in AKT2_4_DE_genes_CRISPRi_positions:
	ax[3,1].axvline(x=gene_pos,color="grey")
ax[3,1].axvline(x = gene_locs.loc[gene_locs["gene"]=="AKT2", 'start'].item(), color = 'blue')
ax[3,1].set_yticks([])
ax[3,1].set_xlabel("Position (bp) along chr19 (hg38)")
ax[3,1].set_title("CRISPRi guide AKT2_4")
ax[3,1].set_xlim(1,chromosome_19_length)

plt.savefig("../figures/AKT2_Dayarmeffectrugplot.png",bbox_inches="tight",dpi=600)




