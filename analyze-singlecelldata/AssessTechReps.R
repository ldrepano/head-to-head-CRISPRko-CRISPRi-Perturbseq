library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(Matrix) 
library(stringr)
library(tidyr)
library(ggsignif)
library(qlcMatrix)
#Assessment of technical replicate correlation:
#Using all cells that pass guide assignment from the indicated sample, 
#cells are pseudobulked by their assigned guide and either their technical 
#replicate or a random label (“pseudoreplicate”). Pearson’s correlation is 
#calculated between the pseudobulked profile of cells carrying the same guide.
#If techreps were a source of noise, correlation of pseudoreps would be higher

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

main<-function(samplename){
  seuratobj<- LoadSeuratRds(paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/",samplename,".Rds",sep=""))
  
  allcells<-Cells(seuratobj)
  seuratobj<-AddMetaData(seuratobj, str_sub(allcells, -1), col.name = "TechRep")
  seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
  var_genes <- VariableFeatures(seuratobj)
  
  
  #pseudobulk all cells that are of the same replicate, carry the same guide 
  pseudobulk<-AggregateExpression(seuratobj, features= var_genes,group.by = c("TechRep","Guide"), return.seurat = TRUE)
  pseudobulked_normalized_data<-pseudobulk@assays$RNA@layers$data
  colnames(pseudobulked_normalized_data)<-Idents(pseudobulk)
  rownames(pseudobulked_normalized_data)<-Features(pseudobulk)
  
  #separate pseudobulk profiles by replicate 
  rep1_colnames<-grep("g1_", colnames(pseudobulked_normalized_data), value = TRUE)
  rep2_colnames<-grep("g2_", colnames(pseudobulked_normalized_data), value = TRUE)
  pseudobulk_rep1 <- pseudobulked_normalized_data[, colnames(pseudobulked_normalized_data) %in% rep1_colnames]
  pseudobulk_rep2 <- pseudobulked_normalized_data[, colnames(pseudobulked_normalized_data) %in% rep2_colnames]
  
  #make sure that columns (guides) of the two replicates are identical 
  colnames(pseudobulk_rep1) <- substring(colnames(pseudobulk_rep1), 4)
  colnames(pseudobulk_rep2) <- substring(colnames(pseudobulk_rep2), 4)
  shared_cols <- intersect(colnames(pseudobulk_rep1), colnames(pseudobulk_rep2))
  pseudobulk_rep1 <- pseudobulk_rep1[, shared_cols]
  pseudobulk_rep2 <- pseudobulk_rep2[, shared_cols]
  
  #correlate the profiles 
  sameguide_crossrep_correlations<-diag(corSparse(pseudobulk_rep1 , pseudobulk_rep2))
  
  #compare sameguide_crossrep_correlations to an arbitrary split within replicates 
  
  seuratobj<-AddMetaData(seuratobj, sample(c("A","B"), length(Cells(seuratobj)),replace=TRUE), col.name = "Pseudoreps")
  
  #pseudobulk all cells that are of the same pseudorep, carry the same guide 
  pseudorep_pseudobulk_seurat<-AggregateExpression(seuratobj, features= var_genes,group.by = c("Pseudoreps","Guide"), return.seurat = TRUE)
  pseudorep_pseudobulk<-pseudorep_pseudobulk_seurat@assays$RNA@layers$data
  colnames(pseudorep_pseudobulk)<-Idents(pseudorep_pseudobulk_seurat)
  rownames(pseudorep_pseudobulk)<-Features(pseudorep_pseudobulk_seurat)
  
  #separate pseudobulk profiles by pseudorep 
  pseudorep1_colnames<-grep("A_", colnames(pseudorep_pseudobulk), value = TRUE)
  pseudorep2_colnames<-grep("B_", colnames(pseudorep_pseudobulk), value = TRUE)
  pseudobulk_pseudorep1 <- pseudorep_pseudobulk[, colnames(pseudorep_pseudobulk) %in% pseudorep1_colnames]
  pseudobulk_pseudorep2 <- pseudorep_pseudobulk[, colnames(pseudorep_pseudobulk) %in% pseudorep2_colnames]
  
  #make sure that columns (guides) of the two replicates are identical 
  colnames(pseudobulk_pseudorep1) <- substring(colnames(pseudobulk_pseudorep1), 3)
  colnames(pseudobulk_pseudorep2) <- substring(colnames(pseudobulk_pseudorep2), 3)
  shared_cols <- intersect(colnames(pseudobulk_pseudorep1), colnames(pseudobulk_pseudorep2))
  pseudobulk_pseudorep1 <- pseudobulk_pseudorep1[, shared_cols]
  pseudobulk_pseudorep2 <- pseudobulk_pseudorep2[, shared_cols]
  
  #correlate the profiles 
  sameguide_crosspseudorep_correlations<-diag(corSparse(pseudobulk_pseudorep1 , pseudobulk_pseudorep2))
  
  
  corr_df<-data.frame(sameguide_crossrep_correlations,sameguide_crosspseudorep_correlations)
  colnames(corr_df)<-c("Technical Replicates","Pseudoreplicates")
  melted_corr_df <- reshape2::melt(corr_df, variable.name = "corr_calculated_across",
                                   value.name = "pearsoncorr")
  
  ggplot(data=melted_corr_df,aes(corr_calculated_across,pearsoncorr))+
    geom_boxplot()+theme_classic()+labs(x="",y="Pearson's Correlation")+
    ggtitle(paste("Guide correlation across technical replicates\n(vs. arbitrary split), ",samplename,sep=""))+
    geom_signif(comparisons = list(c("Technical Replicates", "Pseudoreplicates")),
                test="wilcox.test",map_signif_level = FALSE, textsize = 6)+ylim(0,1)+
    theme(axis.title.y = element_text(size=14),axis.text.x = element_text(size=14))
  ggsave(paste("Figures/TechRepCorrelation/",samplename,".png",sep=""))
}

main("K562_Day7")
main("K562_Day10")

