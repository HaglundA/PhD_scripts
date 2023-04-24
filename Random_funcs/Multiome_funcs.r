

library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(Seurat)
library(Signac)

# read in matrix

#==============================================================================##

# STEP 1: PRE-PROCESSING

#==============================================================================##


# wrapper function to create multiomic seurat object
CreateSeurat_10XMultiome_data=function(h5_mat,
  fragfile,
  min.cells.rna=5,
  min.cells.atac=10){

  message("Reading in data..")
  inputdata<-Read10X_h5(h5_mat)
  rna_counts <- inputdata$`Gene Expression`
  atac_counts <- inputdata$Peaks


  # add in RNA
  message("Initializing Seurat object. Adding in RNA..")
  seuratobj<- CreateSeuratObject(counts = rna_counts,min.cells=min.cells.rna)


  # add in ATAC
  message("Adding in ATAC..")
  options(warn=-1)
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,verbose=F)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"

  frag.file <- fragfile
  chrom_assay <- CreateChromatinAssay(
     counts = atac_counts,
     sep = c(":", "-"),
     genome = 'hg38',
     fragments = frag.file,
     min.cells = min.cells.atac,
     annotation = annotations
   )

  seuratobj[["ATAC"]] <- chrom_assay
  return(seuratobj)
}


# When the Seurat object is created, it creates an element for the fragment file. Because the frag file
# is large, this is just a path and each time you use a function that utilizes the fragment file (for example,
# callPeaks) it expects the path to the fragfile to be correct. If you created the multiomic seurat object on the HPC
# and then move over to the colossus, any function utilising the fragfile won't work. You will have to create a new seurat
# object with the correct frag file path. More info: https://satijalab.org/signac/articles/data_structures.html#the-fragment-class

call_peaks=function(seuratobj,fragpath,macs2.path){
  DefaultAssay(seuratobj)<-"ATAC"
  peaks <- CallPeaks(seuratobj, macs2.path = macs2.path)
  # around 15min of runtime

  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)


  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"

  # quantify counts in each peak: this takes another 15min
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(seuratobj),
    features = peaks,
    cells = colnames(seuratobj)
  )
  #

  # create a new assay using the MACS2 peak set and add it to the Seurat object
  seuratobj[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = fragpath,
    annotation = annotation
  )
  return(seuratobj)
}


DimRed_ATAC=function(seuratobj){

  # this is the OLD Signac tutorial. The NEW signac tutorial suggests calling Peaks
  # with MACS2 and then doing TFIDF on the "Peaks" assay
  DefaultAssay(seuratobj) <- "ATAC"
  seuratobj <- RunTFIDF(seuratobj)
  seuratobj <- FindTopFeatures(seuratobj, min.cutoff = 'q0')
  seuratobj <- RunSVD(seuratobj)
  seuratobj <- RunUMAP(seuratobj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  return(seuratobj)
}


DimRed_RNA=function(seuratobj){
  DefaultAssay(seuratobj) <- "RNA"

  seuratobj <- SCTransform(seuratobj, verbose = FALSE)
  seuratobj<-RunPCA(seuratobj)
  seuratobj=FindNeighbors(seuratobj,dims=1:50)
  seuratobj=FindClusters(seuratobj,pc.use=1:50)
  seuratobj<-RunUMAP(seuratobj,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  return(seuratobj)
}





#

annotate_celltypes=function(seuratobj,cellmarkers){

  wang.markers=readRDS(cellmarkers)
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

  CellTypeTest.wang<-function(x)
  {
    CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05 )]})
    CellType.Names=c(rep("Unclassified",length(HS.FisherTest.zeisel.result)))
    CellType.Names[which(CellType=="Astro")]="Astro"
    CellType.Names[which(CellType=="Endo")]="Endo"
    CellType.Names[which(CellType=="Ex1")]="Excitatory"
    CellType.Names[which(CellType=="Ex2")]="Excitatory"
    CellType.Names[which(CellType=="Ex3e")]="Excitatory"
    CellType.Names[which(CellType=="Ex4")]="Excitatory"
    CellType.Names[which(CellType=="Ex5b")]="Excitatory"
    CellType.Names[which(CellType=="Ex6a")]="Excitatory"
    CellType.Names[which(CellType=="Ex6b")]="Excitatory"
    CellType.Names[which(CellType=="Ex8")]="Excitatory"
    CellType.Names[which(CellType=="Ex9")]="Excitatory"
    CellType.Names[which(CellType=="In1a")]="Inhibitory"
    CellType.Names[which(CellType=="In1b")]="Inhibitory"
    CellType.Names[which(CellType=="In1c")]="Inhibitory"
    CellType.Names[which(CellType=="In3")]="Inhibitory"
    CellType.Names[which(CellType=="In4a")]="Inhibitory"
    CellType.Names[which(CellType=="In4b")]="Inhibitory"
    CellType.Names[which(CellType=="In6a")]="Inhibitory"
    CellType.Names[which(CellType=="In6b")]="Inhibitory"
    CellType.Names[which(CellType=="In7")]="Inhibitory"
    CellType.Names[which(CellType=="In8")]="Inhibitory"
    CellType.Names[which(CellType=="Microglia")]="Microglia"
    CellType.Names[which(CellType=="Oligo")]="Oligo"
    CellType.Names[which(CellType=="OPC")]="OPC"
    CellType.Names[which(CellType=="Per")]="Per"
    CellType.Names
  }

  pbmc=seuratobj
  message("Finding markers per cluster..")
  pbmcmarkers<<-FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threHSold = 0.25,verbose=F)
  # pbmcmarkers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

  HS.MarkersPerCluster.LIST=unstack(pbmcmarkers, pbmcmarkers$gene ~ pbmcmarkers$cluster)
  HS.data.zeisel=pbmc
  HS.uniqueGenesAllClusters=unique(unlist(HS.MarkersPerCluster.LIST))
  HS.FisherTest.zeisel.result=lapply(HS.MarkersPerCluster.LIST,HS.FisherTest.zeisel)
  HS.matrixPval.CellTypeXCluster = matrix(unlist(HS.FisherTest.zeisel.result), ncol = 25, byrow = TRUE)
  colnames(HS.matrixPval.CellTypeXCluster)=names(wang.markers)
  HS.LabelledCellTypeClusters.zeisel=CellTypeTest.wang(HS.matrixPval.CellTypeXCluster)
  names(HS.LabelledCellTypeClusters.zeisel)=names(HS.FisherTest.zeisel.result)
  HS.data.zeisel = RenameIdents(HS.data.zeisel,HS.LabelledCellTypeClusters.zeisel)
  Idents(pbmc)<-Idents(HS.data.zeisel)
  pbmc$RNAannotation<-Idents(pbmc)

  return(pbmc)


}


# This integration function below will perform the integration of ATAC-seq data.
# Note that the FindIntegrationAnchors function uses the "rlsi" reduction parameter:
# this is only available if you use the development version of Seurat. It is the only way
# so far that you can integrate ATAC-seq data.

integrate_samples_ATAC=function(seuratobjlist){

  message("Integrating ATAC first..")
  seuratobjlist<-lapply(seuratobjlist,function(x){
    DefaultAssay(x)<-"ATAC"
    return(x)})

  message("Merging ATAC..")
  mergedobject<-seuratobjlist[[1]]
    for(i in 2:length(seuratobjlist)){
      mergedobject<-merge(mergedobject,seuratobjlist[[i]])
    }

    message("Performing dimred..")
    DefaultAssay(mergedobject)<-"ATAC"
    options(warn=-1)
    mergedobjet<-FindTopFeatures(mergedobject,min.cutoff=10)
    mergedobject <- RunTFIDF(mergedobject)
    mergedobject <- RunSVD(mergedobject)
    mergedobject <- RunUMAP(mergedobject, reduction = "lsi", dims = 2:30)

    mergedobject
    message("finding ATAC integration anchors..")
    integration.anchors <- FindIntegrationAnchors(
    object.list = seuratobjlist,
    anchor.features = rownames(seuratobjlist[[1]]),
    reduction = "rlsi",
    dims = 2:30
  )

     message("Integrating ATAC based on computed embeddings..")

    integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = mergedobject[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:30
  )


  return(integrated)


}


integrate_samples_RNA_sct=function(seuratobjlist){

  seuratobjlist<-lapply(seuratobjlist,SCTransform)
  features=SelectIntegrationFeatures(seuratobjlist,nfeatures=nrow(seuratobjlist[[1]]))

  seuratobjlist<-PrepSCTIntegration(seuratobjlist,anchor.features=features)

  anchors<-FindIntegrationAnchors(seuratobjlist,normalization.method="SCT",anchor.features=features)
  combined<-IntegrateData(anchorset=anchors,normalization.method="SCT")

  combined<-RunPCA(combined,verbose=F)
  combined<-RunUMAP(combined,reduction="pca",dims=1:30)
  return(combined)


}
