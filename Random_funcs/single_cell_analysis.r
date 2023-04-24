

source(single_cell_functions.r)
wang.markers=readRDS("/rds/general/user/ah3918/home/SingleCell/Wang_Markers.rds")


sc_analysis=function(seuratobj,
  normalization.method="LogNormalize",
  scale.factor=10000,
){

  seuratobj<-NormalizeData(seuratobj,
    normalization.method=normalization.method,
    scale.factor=scale.factor)


  seuratobj<-FindVariableFeatures(seuratobj,selection.method="vst")
  seuratobj=ScaleData(seuratobj,features=rownames(seuratobj))

  seuratobj<-RunPCA(seuratobj,npcs=30,verbose=T)

  ElbowPlot(seuratobj)
  DimHeatmap(seuratobj,dims=1:30)

  DimPlot(seuratobj)

  seuratobj<-JackStraw(seuratobj,num.replicate=100)
  seuratobj<-ScoreJackStraw(seuratobj,reduction="pca",dims=1:20)

  JackStrawPlot(seuratobj,dims=1:30)



  seuratobj<-FindNeighbors(seuratobj,dims=30)
  seuratobj<-FindClusters(seuratobj,pc.use=30)
  seuratobj<-RunUMAP(seuratobj,dims=30)
  seuratobj<-RunTSNE(seuratobj,dims=30)

}


# wang.markers=readRDS("/rds/general/user/ah3918/home/SingleCell/Wang_Markers.rds")
#
# #Fishers Exact test for zeisel markers/cluster
#
#
HS.FisherTest.zeisel=function(x)
{
  TMP=matrix(ncol=25,nrow=1)
  HS.uniqueGenesPerCluster=intersect(HS.uniqueGenesAllClusters,x)
  for(j in 1:25)
  {
    TMP.Genes=intersect(HS.uniqueGenesAllClusters,wang.markers[[j]])
    TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,HS.uniqueGenesPerCluster))
    TMP.MAT[1,2]=length(setdiff(TMP.Genes,HS.uniqueGenesPerCluster))
    TMP.MAT[2,1]=length(setdiff(HS.uniqueGenesPerCluster,TMP.Genes))
    TMP.MAT[2,2]=length(HS.uniqueGenesAllClusters)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    TMP[1,j]=fisher.test(TMP.MAT,alternative="greater")$p.value
  }
  TMP
}

# HS.data.zeisel=pbmc
# HS.uniqueGenesAllClusters=unique(unlist(HS.MarkersPerCluster.LIST))
# HS.FisherTest.zeisel.result=lapply(HS.MarkersPerCluster.LIST,HS.FisherTest.zeisel)
# HS.matrixPval.CellTypeXCluster = matrix(unlist(HS.FisherTest.zeisel.result), ncol = 25, byrow = TRUE)
# colnames(HS.matrixPval.CellTypeXCluster)=names(wang.markers)
#
# HS.matrixPval.CellTypeXCluster
#
# CellTypeTest.wang<-function(x)
# {
#   CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05 )]})
#   CellType.Names=c(rep("Unclassified",length(HS.FisherTest.zeisel.result)))
#   CellType.Names[which(CellType=="Astro")]="Astro"
#   CellType.Names[which(CellType=="Endo")]="Endo"
#   CellType.Names[which(CellType=="Ex1")]="Excitatory"
#   CellType.Names[which(CellType=="Ex2")]="Excitatory"
#   CellType.Names[which(CellType=="Ex3e")]="Excitatory"
#   CellType.Names[which(CellType=="Ex4")]="Excitatory"
#   CellType.Names[which(CellType=="Ex5b")]="Excitatory"
#   CellType.Names[which(CellType=="Ex6a")]="Excitatory"
#   CellType.Names[which(CellType=="Ex6b")]="Excitatory"
#   CellType.Names[which(CellType=="Ex8")]="Excitatory"
#   CellType.Names[which(CellType=="Ex9")]="Excitatory"
#   CellType.Names[which(CellType=="In1a")]="Inhibitory"
#   CellType.Names[which(CellType=="In1b")]="Inhibitory"
#   CellType.Names[which(CellType=="In1c")]="Inhibitory"
#   CellType.Names[which(CellType=="In3")]="Inhibitory"
#   CellType.Names[which(CellType=="In4a")]="Inhibitory"
#   CellType.Names[which(CellType=="In4b")]="Inhibitory"
#   CellType.Names[which(CellType=="In6a")]="Inhibitory"
#   CellType.Names[which(CellType=="In6b")]="Inhibitory"
#   CellType.Names[which(CellType=="In7")]="Inhibitory"
#   CellType.Names[which(CellType=="In8")]="Inhibitory"
#   CellType.Names[which(CellType=="Microglia")]="Microglia"
#   CellType.Names[which(CellType=="Oligo")]="Oligo"
#   CellType.Names[which(CellType=="OPC")]="OPC"
#   CellType.Names[which(CellType=="Per")]="Per"
#   CellType.Names
# }
# # Annotate cells based on ATAC only
#
# HS.LabelledCellTypeClusters.zeisel=CellTypeTest.wang(HS.matrixPval.CellTypeXCluster)
# names(HS.LabelledCellTypeClusters.zeisel)=names(HS.FisherTest.zeisel.result)
# HS.data.zeisel = RenameIdents(HS.data.zeisel,HS.LabelledCellTypeClusters.zeisel)
# levels(HS.data.zeisel)
#
# png("11NOV20_UMAP_CellTypes.png")
# options(repr.plot.width=8, repr.plot.height=6)
# DimPlot(HS.data.zeisel, reduction = "umap", label = TRUE, pt.size = 1,label.size=4)
# table(Idents(HS.data.zeisel))
# dev.off()
#
# Idents(pbmc)<-Idents(HS.data.zeisel)
# pbmc$RNAannotation<-Idents(pbmc)
#
# saveRDS(pbmc,"3kHumanBrain_ATACRNA_SeuratObj_annotated.rds")
