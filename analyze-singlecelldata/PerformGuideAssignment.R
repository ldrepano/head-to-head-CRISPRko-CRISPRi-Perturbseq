library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(Matrix) 
library(stringr)

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

main <- function(samplename) {
  guide_reads_per_cell<-read.csv(paste("GuideReadsPerNondyingCell/",samplename,".csv",sep=""),row.names = 1)
  
  guide_purity<-apply(guide_reads_per_cell, 2, function(x) max(x, na.rm = TRUE)) /colSums(guide_reads_per_cell,na.rm=TRUE)
  
  #get data such that each row is a barcode: guide pair 
  guide_reads_per_cell_long<-guide_reads_per_cell
  guide_reads_per_cell_long$GuideID<-rownames(guide_reads_per_cell_long)
  guide_reads_per_cell_long<-reshape2::melt(guide_reads_per_cell_long,id.vars="GuideID",variable.name="CellBarcode",value.name="n_guidereads")
  #smallest # guide reads that could be interpreted as native
  guidereadcutoff<-10
  
  cells_with_guide_purity<- names(guide_purity[guide_purity>0.65])
  maxguidereadcounts<- apply(guide_reads_per_cell, 2, function(x) max(x, na.rm = TRUE))
  cells_with_sufficient_guideguidereads<-names(maxguidereadcounts[maxguidereadcounts>guidereadcutoff])
  datafilt<-guide_reads_per_cell[, c(intersect(cells_with_sufficient_guideguidereads,cells_with_sufficient_guideguidereads))]
  most_abundant_guide <- rownames(datafilt)[apply(datafilt, MARGIN = 2, which.max)]
  guideassignments<- data.frame(Cell=names(datafilt),Guide=most_abundant_guide)
  
  #remove guides with <100 cell coverage 
  cellsperguide<-guideassignments %>% group_by(Guide) %>%tally()
  guides_with_sufficient_coverage<- cellsperguide[cellsperguide$n>=100,]$Guide
  guideassignments_filt<-subset(guideassignments,Guide %in% guides_with_sufficient_coverage)
  write.csv(guideassignments_filt,paste("OneAssignedGuidePerCell_Min100CellsPerGuide/",samplename,sep=""),row.names=FALSE)
  
  cellsperguide<-guideassignments_filt %>% group_by(Guide) %>%tally()
  write.csv(cellsperguide,paste("CellCountsperGuide/",samplename,sep=""),row.names=FALSE)
}

main("K562_Day4")
main("K562_Day7")
main("K562_Day10")
main("K562_Day14_CRISPRko")
main("K562_Day14_CRISPRi")
main("A549_Day7")
