library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(stringr)
setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

#identifies CRISPRko vs CRISPRi enrichment of genes on the target chromosome arm towards telomere 

#process reference info
gene_loc<-read.table("Reference/gene_locations_ensembl115.txt",col.names=c("gene","chr","start","stop"))
gene_to_chr <- setNames( gene_loc$chr,gene_loc$gene)
gene_to_pos <- setNames( gene_loc$start,gene_loc$gene)
#manual correction of targets
gene_to_chr["MOSMO"]<-"16"
gene_to_pos["MOSMO"]<-22007638
gene_to_chr["MGAT4B"]<-"5"
gene_to_pos["MGAT4B"]<-179224597
get_gene_chr<-function(gene){return(gene_to_chr[gene])}

chr_centromere_pos<-read.csv("Reference/UCSC_hg38_chromosome_band_track.csv")
chr_centromere_pos$chr <- gsub("chr", "", chr_centromere_pos$X.chrom)
chr_centromere_pos$arm<-substr(chr_centromere_pos$name,start=1,stop=1)
chr_centromere_pos<-subset(chr_centromere_pos,arm=="p")
chr_centromere_pos<-chr_centromere_pos %>% dplyr::group_by(chr) %>% dplyr::summarise(centromere_pos = max(chromEnd))
chr_to_centromere<- setNames(chr_centromere_pos$centromere_pos,chr_centromere_pos$chr)
get_chr_centromere<-function(chr){return(chr_to_centromere[chr])}

intergenics_info<-read.csv("Reference/intergenic_guides_with_target_site.csv")
#subset to intergenics that target chromosomes that do not harbor any of the target genes for this library
reference_intergenics<-subset(intergenics_info,chromosome %in% c("10","12","22"))
library_info<-read.csv("Reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_short_to_long<-setNames(library_info$guide_id_long, library_info$guide_id_short)
reference_intergenic_guides<-guide_short_to_long[reference_intergenics$guide_name]

get_towardstelomere_enrichment_percell<-function(sample){
  seurat_df<- LoadSeuratRds(paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/K562_Day",sample,".Rds",sep=""))
  seurat_df<-NormalizeData(seurat_df)
  seurat_df<-AddMetaData(seurat_df, sub("_.*", "", seurat_df@meta.data$Guide),col.name="target_gene")
  #remove cells with MPG perturbed (only 3 genes away from centromere)
  seurat_df<-subset(seurat_df,target_gene!="MPG")
  
  target_gene_chromosomes<-unique(lapply(seurat_df@meta.data$target_gene,get_gene_chr))
  #restrict DF to only report data from target gene chromosomes 
  gene_loc<-subset(gene_loc,chr %in% target_gene_chromosomes)
  gene_loc$centromere_pos<-lapply(gene_loc$chr,get_chr_centromere)
  gene_loc$arm<-ifelse(gene_loc$start < gene_loc$centromere_pos , "p", "q")
  gene_to_arm<-setNames(gene_loc$arm, gene_loc$gene)
  gene_to_arm["MOSMO"]<-"p"
  gene_to_arm["MGAT4B"]<-"q"
  
  #only includes genes per day covered by both CRISPRko and CRISPRi
  seurat_CRISPRkosubset_fortargetgenes<-subset(seurat_df,modality=="CRISPRko")
  target_genes_CRISPRko<-unique(seurat_CRISPRkosubset_fortargetgenes@meta.data$target_gene)
  seurat_CRISPRisubset_fortargetgenes<-subset(seurat_df,modality=="CRISPRi")
  target_genes_CRISPRi<-unique(seurat_CRISPRisubset_fortargetgenes@meta.data$target_gene)
  target_genes<-intersect(target_genes_CRISPRko,target_genes_CRISPRi)
  
  downstream_genes<-c()
  for (gene in target_genes){
    gene_chr<-gene_to_chr[gene]
    gene_pos<-gene_to_pos[gene]
    gene_arm<-gene_to_arm[gene]
    downstream_gene_list<-subset(gene_loc,chr==gene_chr)
    downstream_gene_list<-subset(downstream_gene_list,arm==gene_arm)
    if(gene_arm == "p"){
      downstream_gene_list<-subset(downstream_gene_list,start< gene_pos)
    }else{downstream_gene_list<-subset(downstream_gene_list,start> gene_pos)}
    downstream_genes<-c(downstream_genes,list(downstream_gene_list$gene))
  }
  
  enrichment_of_genes_away_from_centromere_CRISPRko<-c()
  enrichment_of_genes_away_from_centromere_CRISPRi<-c()
  enrichment_of_genes_away_from_centromere_REF<-c()
  #essentially tracks # cells per target for each cell set
  target_REF<-c()
  target_CRISPRko<-c()
  target_CRISPRi<-c()
  num_genes_in_enrichment<-c()#number of genes away from centromere that are in seurat object
  # stratify by modality, then compare distribution of enrichment per gene.
  # then repeat for other time points
  for (i in 1:length(target_genes)){
    downstream_genes_subset<-downstream_genes[i]
    seurat_df<-AddModuleScore(seurat_df,features=downstream_genes_subset,name="genes_away_from_centromere")
    num_genes_in_enrichment<-c(num_genes_in_enrichment,length(intersect(downstream_genes[i][[1]],Features(seurat_df))))
    
    reference_intergenics_seurat_obj<-subset(seurat_df,Guide %in% reference_intergenic_guides)
    enrichment_of_genes_away_from_centromere_REF<-c(enrichment_of_genes_away_from_centromere_REF,reference_intergenics_seurat_obj@meta.data$genes_away_from_centromere1)
    target_REF<-c(target_REF,rep(target_genes[i], times = nrow(reference_intergenics_seurat_obj@meta.data)))
    
    seurat_target_subset<-subset(seurat_df,target_gene==target_genes[i])

    CRISPRko_seurat_obj<-subset(seurat_target_subset,modality=="CRISPRko")
    enrichment_of_genes_away_from_centromere_CRISPRko<-c(enrichment_of_genes_away_from_centromere_CRISPRko,CRISPRko_seurat_obj@meta.data$genes_away_from_centromere1)
    target_CRISPRko<-c(target_CRISPRko,CRISPRko_seurat_obj@meta.data$target_gene)

    CRISPRi_seurat_obj<-subset(seurat_target_subset,modality=="CRISPRi")
    enrichment_of_genes_away_from_centromere_CRISPRi<-c(enrichment_of_genes_away_from_centromere_CRISPRi,CRISPRi_seurat_obj@meta.data$genes_away_from_centromere1)
    target_CRISPRi<-c(target_CRISPRi,CRISPRi_seurat_obj@meta.data$target_gene)
    
  }
  print("minimum count of genes towards telomere, in enrichment calculation:")
  print(min(num_genes_in_enrichment))
  
  CRISPRko_df<-as.data.frame(enrichment_of_genes_away_from_centromere_CRISPRko)
  CRISPRko_df<-cbind(CRISPRko_df,as.data.frame(target_CRISPRko))
  colnames(CRISPRko_df) <- c("Enrichment","Target")
  CRISPRko_df$Day<-sample
  CRISPRko_df$modality<-"CRISPRko"
  
  CRISPRi_df<-as.data.frame(enrichment_of_genes_away_from_centromere_CRISPRi)
  CRISPRi_df<-cbind(CRISPRi_df,as.data.frame(target_CRISPRi))
  colnames(CRISPRi_df) <- c("Enrichment","Target")
  CRISPRi_df$Day<-sample
  CRISPRi_df$modality<-"CRISPRi"
  
  REF_df<-as.data.frame(enrichment_of_genes_away_from_centromere_REF)
  REF_df<-cbind(REF_df,as.data.frame(target_REF))
  colnames(REF_df) <- c("Enrichment","Target")
  REF_df$Day<-sample
  REF_df$modality<-"REF"
  
  enrichment_df<- rbind(CRISPRko_df,CRISPRi_df,REF_df)
  
  return(enrichment_df)
}
  
d4_enrichmentdf<-get_towardstelomere_enrichment_percell("4")
d7_enrichmentdf<-get_towardstelomere_enrichment_percell("7")
d10_enrichmentdf<-get_towardstelomere_enrichment_percell("10")
d14_enrichmentdf<-get_towardstelomere_enrichment_percell("14")

#save results since get_towardstelomere_enrichment_percell takes a long time to run
write.csv(d4_enrichmentdf, "TowardsTelomereEnrichment/d4_enrichmentofgenestowardstelomere.csv", row.names = FALSE)
write.csv(d7_enrichmentdf, "TowardsTelomereEnrichment/d7_enrichmentofgenestowardstelomere.csv", row.names = FALSE)
write.csv(d10_enrichmentdf, "TowardsTelomereEnrichment/d10_enrichmentofgenestowardstelomere.csv", row.names = FALSE)
write.csv(d14_enrichmentdf, "TowardsTelomereEnrichment/d14_enrichmentofgenestowardstelomere.csv", row.names = FALSE)

#(reload in if needed)
d4_enrichmentdf<-read.csv("TowardsTelomereEnrichment/d4_enrichmentofgenestowardstelomere.csv")
d7_enrichmentdf<-read.csv("TowardsTelomereEnrichment/d7_enrichmentofgenestowardstelomere.csv")
d10_enrichmentdf<-read.csv("TowardsTelomereEnrichment/d10_enrichmentofgenestowardstelomere.csv")
d14_enrichmentdf<-read.csv("TowardsTelomereEnrichment/d14_enrichmentofgenestowardstelomere.csv")

#zscore enrichment scores to intergenics

get_enrichment_zscore<-function(df){
  
  #per target, get mean&sd enrichment of genes towards telomere in negcon cells
  targets<-unique(df$Target)
  negcon_means<-c()
  negcon_sds<-c()
  negcon_bottomfifthpercentile<-c()
  for(target in targets){
    targetsubset<-subset(df,Target==target)
    targetsubset_negcons<-subset(targetsubset,modality=="REF")
    negcon_enrichment_mean<-mean(targetsubset_negcons$Enrichment)
    negcon_enrichment_sd<-sd(targetsubset_negcons$Enrichment)
    negcon_means<-c(negcon_means,negcon_enrichment_mean)
    negcon_sds<-c(negcon_means,negcon_enrichment_mean)
    
    bottom5thpercentile<-quantile(targetsubset_negcons$Enrichment, probs = 0.05)
    
    negcon_bottomfifthpercentile<-c(negcon_bottomfifthpercentile,
                                              bottom5thpercentile)
  }
  
  target_to_negconmean <- setNames( negcon_means,targets)
  target_to_negconsd <- setNames( negcon_sds,targets)
  target_to_negcon_bottomfifthpercentile <- setNames( negcon_bottomfifthpercentile,targets)
  
  #use negcon means and sds to calculate enrichment zscores
  df$negconmean<-target_to_negconmean[df$Target]
  df$negconsd<-target_to_negconsd[df$Target]
  df$diff_from_negcon_mean<-df$Enrichment-df$negconmean
  df$zscore<-df$diff_from_negcon_mean/df$negconsd
  
  #identify whether or not a cell has lower towards telomere transcript abundance than 95+% negcon cells
  df$negconbottomfifthpercentile<-target_to_negcon_bottomfifthpercentile[df$Target]
  df$is_below_negconbottomfifthpercentile<-df$Enrichment < df$negconbottomfifthpercentile
  
  #clean up columns
  df <- select(df, -c(negconmean, negconsd,diff_from_negcon_mean,Enrichment,negconbottomfifthpercentile))
  
  #remove intergenics
  df<-subset(df,modality!="REF")
  
  return(df)
}

d4_enrichmentdf<-get_enrichment_zscore(d4_enrichmentdf)
d7_enrichmentdf<-get_enrichment_zscore(d7_enrichmentdf)
d10_enrichmentdf<-get_enrichment_zscore(d10_enrichmentdf)
d14_enrichmentdf<-get_enrichment_zscore(d14_enrichmentdf)


AKT2_towardstelomere_genes<-subset(d4_enrichmentdf,Target=="AKT2")

ggplot(data=AKT2_towardstelomere_genes,aes(x=zscore,color=modality))+
  geom_density(linewidth=1.3)+theme_classic()+
  scale_color_manual(values = setNames(c("dodgerblue","limegreen"), c("CRISPRko","CRISPRi")))+
  labs(x="Transcript abundance towards telomere\nZ-scored to negative control cells",
       y="Density",title="Day 4 arm loss in cells with AKT2-targeting guides")
ggsave("Figures/ArmEffects/percent_cells_with_loss_towards_telomere_AKT2_day4.png",width=6,height=3)


get_percent_loss<-function(enrichment_df,sample){
  #get % of towards telomere loss for CRISPRko cells per target
  percent_loss_CRISPRko<-enrichment_df %>%
    dplyr::group_by(Target, modality) %>%
    dplyr::summarise(
      percent_below = 100*mean(is_below_negconbottomfifthpercentile))%>%
    subset(modality=="CRISPRko")
  percent_loss_CRISPRko$Day<-sample
  
  #repeat for CRISPRi
  percent_loss_CRISPRi<-enrichment_df %>%
    dplyr::group_by(Target, modality) %>%
    dplyr::summarise(
      percent_below = 100*mean(is_below_negconbottomfifthpercentile))%>%
    subset(modality=="CRISPRi")
  percent_loss_CRISPRi$Day<-sample
  
  return(list(percent_loss_CRISPRko,percent_loss_CRISPRi))
}


d4_lossrate<-get_percent_loss(d4_enrichmentdf,sample="4")
d4_lossrate_CRISPRko<-d4_lossrate[[1]]
d4_lossrate_CRISPRi<-d4_lossrate[[2]]
d7_lossrate<-get_percent_loss(d7_enrichmentdf,sample="7")
d7_lossrate_CRISPRko<-d7_lossrate[[1]]
d7_lossrate_CRISPRi<-d7_lossrate[[2]]
d10_lossrate<-get_percent_loss(d10_enrichmentdf,sample="10")
d10_lossrate_CRISPRko<-d10_lossrate[[1]]
d10_lossrate_CRISPRi<-d10_lossrate[[2]]
d14_lossrate<-get_percent_loss(d14_enrichmentdf,sample="14")
d14_lossrate_CRISPRko<-d14_lossrate[[1]]
d14_lossrate_CRISPRi<-d14_lossrate[[2]]

#here I am distinguishing target genes by color and shape 
cats <- unique(d4_lossrate_CRISPRko$Target)
n_colors<-12
colors <- scales::hue_pal()(n_colors)
pretty_palette_12 <- c(
  "#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
  "#999999","#8DD3C7","#BEBADA", "#FB8072","#80B1D3")
shapes <- 16:18 
df_map <- data.frame(Target = cats) %>%
  dplyr::mutate(
    idx = row_number() - 1,
    color_id = (idx %% n_colors) + 1,
    shape_id = (idx %/% n_colors) + 1,
    color = pretty_palette_12[color_id],
    shape = shapes[shape_id])

alldays_lossrate_CRISPRko<-do.call(rbind,list(d4_lossrate_CRISPRko,d7_lossrate_CRISPRko,d10_lossrate_CRISPRko,d14_lossrate_CRISPRko))
alldays_lossrate_CRISPRi<-do.call(rbind,list(d4_lossrate_CRISPRi,d7_lossrate_CRISPRi,d10_lossrate_CRISPRi,d14_lossrate_CRISPRi))
alldays_lossrate<-rbind(alldays_lossrate_CRISPRko,alldays_lossrate_CRISPRi)
alldays_lossrate$Day <- factor(alldays_lossrate$Day, levels = c("4", "7", "10","14"))
ggplot(alldays_lossrate,aes(Day, percent_below,dodge=modality))+
  geom_boxplot(aes(Day, percent_below,linetype=modality),outlier.shape = NA)+
  theme_classic()+
  #geom_jitter(aes(Day, percent_below,color=Target,shape=Target),width=0.1)+
  geom_point(aes(Day, percent_below,dodge=modality,color=Target,shape=Target),position = position_jitterdodge(jitter.width = 0.2))+
  scale_y_continuous(breaks = seq(0, 45, by = 5))+
  scale_color_manual(values = setNames(df_map$color, df_map$Target)) +
  scale_shape_manual(values = setNames(df_map$shape, df_map$Target))+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 14))+
  labs(title="% cells per target with loss\ntowards telomere on target chromosome arm",
       y= "% cells per target")+
  stat_compare_means(label="p.format",method="wilcox.test")+
  coord_cartesian(ylim = c(0,28))
ggsave("Figures/ArmEffects/percent_cells_with_loss_towards_telomere.png",width=8,height=5)

mean(d4_lossrate_CRISPRko$percent_below) 

