library(Seurat)
library(stringr)
setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

# get # cells after removing those with >15% mitochondrial transcripts
get_assignment_rate<-function(samplename){
  print(samplename)
  seuratobj_postqc<-LoadSeuratRds(paste("PreassignmentReplicatesCombinedNondyingCellsRds/",samplename,".Rds",sep=""))
  num_passqc_cells<-length(Cells(seuratobj_postqc))
  print("# cells that pass DRAGEN qc and have <=15% mitochondrial transcripts")
  print(num_passqc_cells)
  guideassignments<- read.csv(paste("CellCountsperGuide/",samplename,sep=""))
  num_assigned_cells<-sum(guideassignments$n)
  print("# cells remaining after guide assignment, removing targets with <100 cells per guide")
  print(num_assigned_cells)
  print("assignment rate:")
  print(num_assigned_cells/num_passqc_cells)
}

get_assignment_rate("K562_Day4")
get_assignment_rate("K562_Day7")
get_assignment_rate("K562_Day10")
get_assignment_rate("K562_Day14_CRISPRko")
get_assignment_rate("K562_Day14_CRISPRi")
get_assignment_rate("A549_Day7")


#----------------get expression of library targets in intergenics --------------------

get_intergenics<-function(samplename){
  seuratobj_postqc<-LoadSeuratRds(paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/",samplename,".Rds",sep=""))
  seuratobj_subset_to_intergenics<-subset(seuratobj_postqc,grepl("INTERGENIC", Guide))
  seuratobj_subset_to_intergenics<-NormalizeData(seuratobj_subset_to_intergenics)
}

intergenics_K562Day4<-get_intergenics("K562_Day4")
intergenics_K562Day7<-get_intergenics("K562_Day7")
intergenics_K562Day10<-get_intergenics("K562_Day10")
intergenics_K562Day14_CRISPRko<-get_intergenics("K562_Day14_CRISPRko")
intergenics_K562Day14_CRISPRi<-get_intergenics("K562_Day14_CRISPRi")
intergenics_A549Day7<-get_intergenics("A549_Day7")

allsamples<-merge(x = intergenics_K562Day4,y = c(intergenics_K562Day7,
                                                 intergenics_K562Day10,intergenics_K562Day14_CRISPRko,
                                                 intergenics_K562Day14_CRISPRi,intergenics_A549Day7),
                  add.cell.ids<-c("Day4","Day7","Day10","Day14 CRISPRko","Day14 CRISPRi", "Day7 A549"),
                  merge.data=TRUE)


allsamples$sample<- str_split_i(Cells(allsamples), "_", 1)
Idents(allsamples)<-"sample"
Idents(allsamples) <- rev(levels(Idents(allsamples)))

#sort library targets by average expression across samples
guide_membership<- read.csv("Reference/CRISPRko-CRISPRi-perturbseq-benchmark-guides.csv")
guide_membership$target<-sub("_.*", "", guide_membership$guide_id_short)
library_targets<-unique(guide_membership$target)
avg.exp <- AverageExpression(object = allsamples, features=library_targets)
avg.exp.matrix <- avg.exp$RNA  
gene.means <- rowMeans(avg.exp.matrix)
library_targets_sorted <- sort(gene.means, decreasing = TRUE)

DotPlot(allsamples,features=names(library_targets_sorted),
        group.by="orig.ident",
        scale=FALSE)+
  labs(y="",x="",title="Library targets in negative control cells")+RotatedAxis()+
  scale_color_distiller(palette="GnBu",direction = 0)
ggsave("Figures/intergenics_librarytargets_expression_dotplot.png",height=6,width=10)

#------------UMAP------------

get_UMAPs<-function(samplename){
  seuratobj_postqc<-LoadSeuratRds(paste("NondyingCellsWholeTranscriptomeWithGuideAssignments/",samplename,".Rds",sep=""))
  seuratobj_postqc<-AddMetaData(seuratobj_postqc, sub("_.*", "", seuratobj_postqc@meta.data$Guide),col.name="target_gene")
  seuratobj_postqc<-NormalizeData(seuratobj_postqc)
  seuratobj_postqc<-ScaleData(seuratobj_postqc)
  seuratobj_postqc<-FindVariableFeatures(seuratobj_postqc)
  seuratobj_postqc <- RunPCA(seuratobj_postqc)
  seuratobj_postqc <- RunUMAP(seuratobj_postqc, dims = 1:10, n.neighbors=8)
  seuratobj_postqc$is_CRISPRko<-seuratobj_postqc$modality=="CRISPRko"
  seuratobj_postqc$is_CRISPRko <- factor(seuratobj_postqc$is_CRISPRko, levels = c(TRUE,FALSE))
  Idents(object = seuratobj_postqc) <- 'is_CRISPRko'
  CRISPRkoumap<-DimPlot(seuratobj_postqc, reduction = "umap",cols = c("navy", "grey"))+labs(color = "Cell perturbed\nwith CRISPRko")+
    ggtitle(samplename)
  seuratobj_postqc$is_CRISPRi<-seuratobj_postqc$modality=="CRISPRi"
  seuratobj_postqc$is_CRISPRi <- factor(seuratobj_postqc$is_CRISPRi, levels = c(TRUE,FALSE))
  Idents(object = seuratobj_postqc) <- 'is_CRISPRi'
  CRISPRiumap<-DimPlot(seuratobj_postqc, reduction = "umap",cols = c("chartreuse2", "grey"))+labs(color = "Cell perturbed\nwith CRISPRi")+
    ggtitle(samplename)
  Idents(object = seuratobj_postqc) <- "target_gene"
  targetumap<-DimPlot(seuratobj_postqc, reduction = "umap")+labs(color = "Target of guide in cell")+
    ggtitle(samplename)
  return(list(CRISPRkoumap, CRISPRiumap,targetumap))
}

d4_UMAPs<-get_UMAPs("K562_Day4")
d4_UMAPs[1][[1]]+d4_UMAPs[2][[1]]+d4_UMAPs[3][[1]]
ggsave("Figures/UMAPs/Day4_UMAPs.png",width=14,height=4)

d7_UMAPs<-get_UMAPs("K562_Day7")
d7_UMAPs[1][[1]]+d7_UMAPs[2][[1]]+d7_UMAPs[3][[1]]
ggsave("Figures/UMAPs/Day7_UMAPs.png",width=14,height=4)

d10_UMAPs<-get_UMAPs("K562_Day10")
d10_UMAPs[1][[1]]+d10_UMAPs[2][[1]]+d10_UMAPs[3][[1]]
ggsave("Figures/UMAPs/Day10_UMAPs.png",width=14,height=4)

d14_UMAPs_CRISPRko<-get_UMAPs("K562_Day14_CRISPRko")
d14_UMAPs_CRISPRi<-get_UMAPs("K562_Day14_CRISPRi")
d14_UMAPs_CRISPRi[3][[1]]+d14_UMAPs_CRISPRko[3][[1]]
ggsave("Figures/UMAPs/Day14_UMAPs.png",width=14,height=4)



