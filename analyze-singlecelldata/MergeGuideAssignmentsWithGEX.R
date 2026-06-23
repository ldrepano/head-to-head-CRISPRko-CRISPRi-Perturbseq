library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(Matrix) 
library(stringr)

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

guide_membership<- read.csv("Reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
CRISPRi_guides<-subset(guide_membership,modality=="CRISPRi")$guide_id_long

main <- function(samplename) {
  seuratobj<- LoadSeuratRds(paste("PreassignmentReplicatesCombinedNondyingCellsRds/",samplename,".Rds",sep=""))
  guideassignments<- read.csv(paste("OneAssignedGuidePerCell_Min100CellsPerGuide/",samplename,sep=""))
  rownames(guideassignments)<-guideassignments$Cell
  guideassignments$Cell<- NULL
  seuratobj <- AddMetaData(object = seuratobj,metadata = guideassignments)
  seuratobj_assigned<-subset(seuratobj, is.na(Guide)==FALSE) 
  
  seuratobj_assigned$modality <- ifelse(seuratobj_assigned$Guide %in% CRISPRi_guides, "CRISPRi", "CRISPRko")
  seuratobj_assigned$modality <- ifelse(str_detect(seuratobj_assigned$Guide,"CONTROL"), "Both", seuratobj_assigned$modality)
  
  SaveSeuratRds(seuratobj_assigned,file = paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/",samplename,".Rds",sep=""))
}

main("K562_Day4")
main("K562_Day7")
main("K562_Day10")
main("K562_Day14_CRISPRko")
main("K562_Day14_CRISPRi")
main("A549_Day7")

# merge Day 14 CRISPRko and CRISPRi data 
d14_K562_CRISPRko<- LoadSeuratRds(paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/K562_Day14_CRISPRko.Rds",sep=""))
d14_K562_CRISPRi<- LoadSeuratRds(paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/K562_Day14_CRISPRi.Rds",sep=""))
both_d14<-merge(d14_K562_CRISPRko,d14_K562_CRISPRi)
SaveSeuratRds(both_d14,file = paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/K562_Day14.Rds",sep=""))
