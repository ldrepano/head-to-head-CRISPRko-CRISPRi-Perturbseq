# using Mixscale to analyze perturbation effects, the degree of guide agreement
# just showing day7 as proof of concept 

#installations I had to perform to run Mixscale
#devtools::install_github('immunogenomics/presto')
#install.packages("PMA")
#install.packages("protoclust")
#BiocManager::install("glmGamPoi")
#devtools::install_github("satijalab/Mixscale")

library(Mixscale)
library(Seurat)

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")
d7_seurat<- LoadSeuratRds("NondyingCellsWholeTranscriptomeWithGuideAssignments/K562_Day7.Rds")


guide_membership<- read.csv("Reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guidename_dict <- setNames(guide_membership$guide_id_short, guide_membership$guide_id_long)
short_guide_names<-guidename_dict[d7_seurat@meta.data$Guide]
names(short_guide_names)<-NULL
d7_seurat<-AddMetaData(d7_seurat,short_guide_names,col.name="guideID_short")

CRISPRi<- subset(d7_seurat,modality %in% c("CRISPRi","Both"))
CRISPRko<- subset(d7_seurat,modality %in% c("CRISPRko","Both"))

get_mixscaleresults_forplot<-function(sample){
  
  sample<-NormalizeData(sample)
  sample<-ScaleData(sample)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  sample<-RunPCA(sample)
  
  sample<-AddMetaData(sample, sub("_.*", "", sample@meta.data$Guide),col.name="target_gene")
  
  #identified distance from intergenic nearest neighbors 
  with_mixscale <- CalcPerturbSig(
    object = sample, 
    assay = "RNA", 
    slot = "data", 
    gd.class ="target_gene", 
    nt.cell.class = "INTERGENIC", 
    reduction = "pca", 
    ndims = 10, 
    num.neighbors = 20, 
    split.by = NULL) 
  
  #gets perturbation score for each cell
  with_mixscale <-RunMixscale(
    object = with_mixscale, 
    assay = "PRTB", 
    slot = "scale.data", 
    labels = "target_gene", 
    nt.class.name = "INTERGENIC", 
    min.de.genes = 5, 
    logfc.threshold = 0,
    de.assay = "RNA",
    max.de.genes = 100, 
    new.class.name = "mixscale_score", 
    fine.mode = F, 
    verbose = F, 
    split.by = NULL)
  
  
  #subset data to exclude the intergenics from the plot ()
  return(with_mixscale)
}

CRISPRi_with_mixscale<-get_mixscaleresults_forplot(CRISPRi)
SaveSeuratRds(CRISPRi_with_mixscale,file = "MixscaleResults/K562_Day7_CRISPRi.Rds")
CRISPRko_with_mixscale<-get_mixscaleresults_forplot(CRISPRko)
SaveSeuratRds(CRISPRko_with_mixscale,file = "MixscaleResults/K562_Day7_CRISPRko.Rds")

#repeat for A549
A549_d7_seurat<- LoadSeuratRds("NondyingCellsWholeTranscriptomeWithGuideAssignments/A549_Day7.Rds")
short_guide_names<-guidename_dict[A549_d7_seurat@meta.data$Guide]
names(short_guide_names)<-NULL
A549_d7_seurat<-AddMetaData(A549_d7_seurat,short_guide_names,col.name="guideID_short")
A549_CRISPRi<- subset(A549_d7_seurat,modality %in% c("CRISPRi","Both"))
A549_CRISPRi_with_mixscale<-get_mixscaleresults_forplot(A549_CRISPRi)
SaveSeuratRds(A549_CRISPRi_with_mixscale,file = "MixscaleResults/A549_Day7_CRISPRi.Rds")
