library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(Matrix) 
library(stringr)
library(tidyr)

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

main <- function(samplename) {
  seuratobj<- LoadSeuratRds(paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/",samplename,".Rds",sep=""))
  
  # getting grna_target_data_frame with guide, target column  
  assignments<- read.csv(paste("OneAssignedGuidePerCell_Min100CellsPerGuide/",samplename,sep=""))
  grna_target_data_frame<- assignments
  grna_target_data_frame$grna_id<-grna_target_data_frame$Guide
  #call all Intergenics non-targeting so that SCEPTRE recognizes them as negative controls
  intergenic_control_guides<-unique(grna_target_data_frame$grna_id[grepl("CONTROL",grna_target_data_frame$grna_id)])
  grna_target_data_frame$grna_target<-grna_target_data_frame$grna_id
  grna_target_data_frame[grna_target_data_frame$grna_id %in% intergenic_control_guides, ]$grna_target<-"non-targeting"
  grna_target_data_frame$grna_target<- str_replace_all(grna_target_data_frame$grna_target, "_[A-Z]*", "")
  grna_target_data_frame$Cell<- NULL
  grna_target_data_frame$Guide<- NULL
  grna_target_data_frame<- grna_target_data_frame[!duplicated(grna_target_data_frame), ]
  
  #identifying guide assignments for all cells in experiment 
  #forcing SCEPTRE to use determined assignments (since low moi defaults to highest count assignment approach)
  assignments$val<-1
  #format: guide:cell pair is 1 if it is the assignment
  assignments_wide<- as.data.frame(pivot_wider(assignments,names_from = Cell, values_from = val))
  # guide:cell pairs for which it is not the assignment have a value of 0 (regardless of actual UMI count)
  assignments_wide<-replace(assignments_wide, is.na(assignments_wide), 0)
  rownames(assignments_wide)<-assignments_wide$Guide
  assignments_wide$Guide<- NULL
  grna_matrix<- as.matrix(assignments_wide)
  
  
  response_df<- as.data.frame(seuratobj@assays$RNA@layers$counts)
  colnames(response_df)<- rownames(data.frame(seuratobj@active.ident))
  rownames(response_df)<-rownames(as.data.frame(seuratobj@assays$RNA@features))
  #note: supplying column names so that they don't need to be in the same order as the guide matrix
  response_matrix<- as.matrix(response_df)
  rm(response_df)
  rm(seuratobj)
  gc()
  #edit response matrix and guide matrix to only contain cells in both
  assigned_guide_cells<- intersect(colnames(response_matrix),colnames(grna_matrix))
  response_matrix<- response_matrix[,assigned_guide_cells]
  grna_matrix<- grna_matrix[,assigned_guide_cells]

  
  
  write.table(grna_target_data_frame,paste("SCEPTREInputs/",samplename,"grna_target_data_frame.txt",sep=""),row.names = FALSE,quote=FALSE)
  rm(grna_target_data_frame)
  writegz <- gzfile(paste("SCEPTREInputs/",samplename,"response_matrix.txt.gz",sep=""), open = "w")
  write.table(response_matrix, writegz,quote=FALSE)
  close(writegz)
  rm(response_matrix)
  gc()
  write.table(grna_matrix,paste("SCEPTREInputs/",samplename,"grna_matrix.txt",sep=""),quote=FALSE)
  rm(grna_matrix)
}

main("K562_Day4")
main("K562_Day7")
main("K562_Day10")
main("K562_Day14_CRISPRko")
main("K562_Day14_CRISPRi")
main("A549_Day7")
