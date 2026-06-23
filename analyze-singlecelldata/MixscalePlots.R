library(Mixscale)
library(Seurat)

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

CRISPRi_with_mixscale<-LoadSeuratRds("MixscaleResults/K562_Day7_CRISPRi.Rds")
CRISPRko_with_mixscale<-LoadSeuratRds("MixscaleResults/K562_Day7_CRISPRko.Rds")

CRISPRi_with_mixscale_STK24<- subset(CRISPRi_with_mixscale, subset = target_gene == "STK24")
RidgePlot(
  CRISPRi_with_mixscale_STK24,
  features = "mixscale_score", 
  group.by = "guideID_short") + NoLegend()+labs(title="CRISPRi Day 7",x="Mixscale Perturbation Score",y="Guide")
ggsave("Figures/Mixscale/day7_CRISPRi_perturbation_score_distribution_STK24.png", width = 7, height = 4)  

CRISPRko_with_mixscale_STK24<- subset(CRISPRko_with_mixscale, subset = target_gene == "STK24")
RidgePlot(
  CRISPRko_with_mixscale_STK24,
  features = "mixscale_score", 
  group.by = "guideID_short") + NoLegend()+labs(title="CRISPRko Day 7",x="Mixscale Perturbation Score",y="Guide")
ggsave("Figures/Mixscale/day7_CRISPRko_perturbation_score_distribution_STK24.png", width = 7, height = 4)  

CRISPRi_with_mixscale_CTBP1<- subset(CRISPRi_with_mixscale, subset = target_gene == "CTBP1")
RidgePlot(
  CRISPRi_with_mixscale_CTBP1,
  features = "mixscale_score", 
  group.by = "guideID_short") + NoLegend()+labs(title="CRISPRi Day 7",x="Mixscale Perturbation Score",y="Guide")
ggsave("Figures/Mixscale/day7_CRISPRi_perturbation_score_distribution_CTBP1.png", width = 7, height = 4)  

#repeat CRISPRi STK24 for A459
CRISPRi_with_mixscale<-LoadSeuratRds("MixscaleResults/A549_Day7_CRISPRi.Rds")
CRISPRi_with_mixscale_STK24<- subset(CRISPRi_with_mixscale, subset = target_gene == "STK24")
RidgePlot(
  CRISPRi_with_mixscale_STK24,
  features = "mixscale_score", 
  group.by = "guideID_short") + NoLegend()+labs(title="CRISPRi Day 7 A549",x="Mixscale Perturbation Score",y="Guide")
ggsave("Figures/Mixscale/A549_day7_CRISPRi_perturbation_score_distribution_STK24.png", width = 7, height = 4)  



