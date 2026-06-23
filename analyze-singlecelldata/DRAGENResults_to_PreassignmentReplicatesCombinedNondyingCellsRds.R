library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(Matrix) 
library(stringr)

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

main <- function(samplename) {
  #no technical replicates at Day 14, and Day 4 rep1 removed for low depth
  if (str_detect(samplename, "Day4") || str_detect(samplename, "Day14")){
    data <- Read10X(data.dir = paste("DRAGENResults/",samplename,sep=""))
    seuratobj <- CreateSeuratObject(counts = data$`Gene Expression`, project = samplename, min.cells = 3)
  } else {
    Rep1data <- Read10X(data.dir = paste("DRAGENResults/",samplename,"_Rep1",sep=""))
    seuratobj_Rep1 <- CreateSeuratObject(counts = Rep1data$`Gene Expression`, project = samplename, min.cells = 3)
    
    Rep2data <- Read10X(data.dir = paste("DRAGENResults/",samplename,"_Rep2",sep=""))
    seuratobj_Rep2 <- CreateSeuratObject(counts = Rep2data$`Gene Expression`, project = samplename, min.cells = 3)
    
    seuratobj <- merge(seuratobj_Rep1, y = seuratobj_Rep2,project = samplename)
    seuratobj[["RNA"]] <- JoinLayers(seuratobj[["RNA"]])
  }
  
  gc() #opportunity to free unused memory
  
  #get % mitochondrial reads in samples
  seuratobj[["percent.mt"]] <- PercentageFeatureSet(seuratobj, pattern = "^MT-")
  
  #retain only cells for which <15% of transcripts are mitochondrial
  MT_highcutoff<- 15 
  
  qc_plots<- ggplot(mapping = aes(seuratobj$nFeature_RNA)) + 
    geom_histogram(fill="skyblue",binwidth=100) + theme_classic()+
    labs(x = "# Genes per Cell") + 
    ggplot(mapping = aes(seuratobj$nCount_RNA)) + 
    geom_histogram(fill="skyblue",binwidth=200) + theme_classic()+
    labs(x = "# UMIs per Cell",title=samplename)+
    theme(plot.title = element_text(hjust = 0.5))+
    ggplot(mapping = aes(seuratobj$percent.mt)) + 
    geom_histogram(fill="skyblue",binwidth=0.5) + theme_classic()+
    labs(x = "% mitochondrial transcripts")+
    geom_vline(xintercept=MT_highcutoff, linetype="dashed", color="red")+
    xlim(0,20)
  ggsave(paste("Figures/SeuratQC/",samplename,".png",sep=""), width = 15, height = 4)
  #filter data 
  seuratobj_filt<- subset(seuratobj,percent.mt < MT_highcutoff )
  
  SaveSeuratRds(seuratobj_filt,file = paste("PreassignmentReplicatesCombinedNondyingCellsRds/",samplename,".Rds",sep=""))
  
  #get guide info filtered to barcode
  if (str_detect(samplename, "Day4") || str_detect(samplename, "Day14")){
    guides<- data$`CRISPR Direct Capture`
    guides_df<-data.frame(guides)
  } else{
    guides_Rep1<- Rep1data$`CRISPR Direct Capture`
    #add suffix to match cell barcode from seurat merge
    colnames(guides_Rep1) <- paste(colnames(guides_Rep1), "1", sep = "_")
    guides_df_Rep1<-data.frame(guides_Rep1)
    guides_Rep2<- Rep2data$`CRISPR Direct Capture`
    colnames(guides_Rep2) <- paste(colnames(guides_Rep2), "2", sep = "_")
    guides_df_Rep2<-data.frame(guides_Rep2)
    guides_df <- cbind(guides_df_Rep1, guides_df_Rep2)
  }
  
  postfiltering_barcodes<- rownames(seuratobj_filt@assays[["RNA"]]@cells)
  guides_df_filt<- guides_df[postfiltering_barcodes]
  guide_UMIs_per_cell<- colSums(guides_df_filt)
  write.csv(guides_df_filt,paste("GuideReadsPerNondyingCell/",samplename,".csv",sep=""))
}

main("K562_Day4")
main("K562_Day7")
main("K562_Day10")
main("K562_Day14_CRISPRko")
main("K562_Day14_CRISPRi")
main("A549_Day7")



