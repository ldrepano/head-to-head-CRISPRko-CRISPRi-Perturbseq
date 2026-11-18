#11/17/26
#plotting results of GSEA performed with KEGG pathway gene sets on DESeq2-computed differential expression of CRISPRko-perturbed cells relative
    #to CRISPRi perturbed cells

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.rc('pdf', fonttype=42)

#read in gsea results
gsea = pd.read_table('GSEA_results/DESeq2_GSEA_Day7_CellLinesMerged.txt')

#define gene sets of interest
ddr_sets = ['Mismatch repair','Homologous recombination','Non-homologous end-joining','Base excision repair','Nucleotide excision repair',
           'Fanconi anemia pathway','p53 signaling pathway','Cell cycle']

#subset to gene sets of interest
gsea_ddr = gsea[gsea['Description'].isin(ddr_sets)]


#barplot + star next to the significant ones
plot_df = gsea_ddr.melt(id_vars='Description',value_vars=['NES_K562','NES_A549','FDR_K562','FDR_A549'],var_name='metric')
plot_df[['metric','Cell line']] = plot_df['metric'].str.split('_',expand=True)
plot_df['Cell line']=plot_df['Cell line'].replace({'A549':'A549 (p53 wild-type)','K562':'K562 (p53 null)'})
plot_df = plot_df.pivot(index=['Description','Cell line'],columns='metric',values='value').reset_index()

#sorting such that the gene sets appear in order of NES score
plot_df = plot_df.sort_values(by='Cell line').sort_values(by='NES')
sns.barplot(data=plot_df,y='Description',x='NES',hue='Cell line',
            palette={'A549 (p53 wild-type)':'blue','K562 (p53 null)':'gray'})

plot_df['star_xpos']=np.where(plot_df['NES']<0,plot_df['NES']-0.2,plot_df['NES']+0.2)
plot_df['star_xpos']=np.where(plot_df['FDR']>=0.1,np.nan,plot_df['star_xpos'])
sns.stripplot(data=plot_df,x='star_xpos',y='Description',marker='*',size=15,hue='Cell line',dodge=True,legend=False,jitter=False,
              palette={'A549 (p53 wild-type)':'blue','K562 (p53 null)':'gray'})

plt.legend(bbox_to_anchor=(1,1))
plt.title('GSEA (with KEGG gene sets) of DDR pathways in CRISPRko vs CRISPRi cells, Day 7')
plt.ylabel('')
plt.xlabel('Normalized Enrichment Score')

plt.savefig("../figures/gsea_results_ddr_pathways_ko_vs_i_with_deseq2_day7_K562_vs_A549_barplot.png",bbox_inches="tight",dpi=600)
plt.savefig("../figures/gsea_results_ddr_pathways_ko_vs_i_with_deseq2_day7_K562_vs_A549_barplot.pdf",bbox_inches="tight",dpi=600)