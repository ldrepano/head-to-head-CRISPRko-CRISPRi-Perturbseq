library(Seurat)
library(stringr)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(msigdbr)
library("org.Hs.eg.db", character.only = TRUE)
setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

get_de_results<-function(samplename){
  seuratobj<-LoadSeuratRds(paste0('NondyingCellsWholeTranscriptomeWithGuideAssignments/',samplename,'.Rds'))
  seuratobj<-subset(seuratobj,modality=='CRISPRko'|modality=='CRISPRi')
  seuratobj<-AddMetaData(seuratobj, sub("_.*", "", seuratobj$Guide),col.name="target_gene")
  pseudobulk<-AggregateExpression(seuratobj,
                                  assays='RNA',
                                  group.by = c("target_gene","modality"),
                                  return.seurat = FALSE)$RNA
  pseudobulk<-pseudobulk[rowSums(pseudobulk>=10)>10,]
  
  #format coldata
  coldata <- data.frame(sample = colnames(pseudobulk))
  coldata$target_gene<-str_split_i(coldata$sample, "_", 1)
  coldata$modality<-str_split_i(coldata$sample, "_", 2)
  
  #run deseq2
  dds <- DESeqDataSetFromMatrix(countData = pseudobulk,
                                colData = coldata,
                                design = ~ target_gene + modality)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("modality","CRISPRko","CRISPRi"))
  res<-as.data.frame(res)
  res<-res%>%arrange(desc(log2FoldChange))
  return(res)
}

de_results_k562<-get_de_results('K562_Day7')
de_results_a549<-get_de_results('A549_Day7')

get_gsea_results_kegg<-function(de_results){
  converted_ids <- bitr(rownames(de_results), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  converted_ids <- converted_ids %>% distinct(ENTREZID, .keep_all = TRUE)
  de_results$SYMBOL<-rownames(de_results)
  de_results_with_converted_ids<-inner_join(de_results,converted_ids,by = "SYMBOL")
  gsea_input<-setNames(de_results_with_converted_ids$log2FoldChange, de_results_with_converted_ids$ENTREZID)
  gsea_input<-na.omit(gsea_input)
  
  gse<- gseKEGG(geneList     = gsea_input,
                organism     = "hsa",
                minGSSize    = 3,
                maxGSSize    = 800,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                keyType       = "ncbi-geneid",
                seed = TRUE)
  return(gse)
}

gse_k562<-get_gsea_results_kegg(de_results_k562)
gse_a549<-get_gsea_results_kegg(de_results_a549)

gse_k562_res<-dplyr::rename(gse_k562@result,NES_K562=NES,FDR_K562=p.adjust)
gse_a549_res<-dplyr::rename(gse_a549@result,NES_A549=NES,FDR_A549=p.adjust)

gse_res_merged<-inner_join(gse_k562_res,gse_a549_res,by=c('ID','Description'))
write.table(gse_res_merged,file="CRISPRko_vs_CRISPRi_global/DESeq2_GSEA_Day7_CellLinesMerged.txt",
            sep = "\t", row.names = FALSE)

