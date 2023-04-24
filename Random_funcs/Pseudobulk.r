

# This function assumes:

# 1. Input are dgcMatrix.rds files
# 2. Cell types are labelled
# 3. Covariates have been added in the meta data
# 4. You have done all pre-processing steps (doubletfinder,emptydrops, MT filter)

library(Seurat)

pseudobulk_scex=function(inputfiles){


  seuratlist<-lapply(inputfiles,readRDS)

  mergeseurat<-seuratlist[[1]]

  # SCTransform is done on the merged Seurat object
  message("merging seurat objects..")
  for(i in 2:length(seuratlist)){
    mergeseurat<-merge(mergeseurat,seuratlist[[i]])
  }

  options(warn=-1)
  message("Running SCTransform..")
  mergeseurat<-SCTransform(mergeseurat,variable.features.n=nrow(mergeseurat),vars.to.regress=c("Age","PMI","Gender"),verbose=FALSE)
  splitlist<-SplitObject(mergeseurat,split.by="orig.ident")
  samplenames<-names(splitlist)


  celltype_means_SCTransformed=function(seuratobj){
    # gets the mean per gene, per cell type, per individual (although this is done on
    # the Pearson Residuals from the SCTransformed data)
    library(Seurat)
    tmp=list(oligos=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Oligo")),assay='SCT',slot='scale.data'),
    astros=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Astro")),assay='SCT',slot='scale.data'),
    Exc_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Excitatory")),assay='SCT',slot='scale.data'),
    Inh_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Inhibitory")),assay='SCT',slot='scale.data'),
    microglia=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Microglia")),assay='SCT',slot='scale.data'))

    tmp=lapply(tmp,rowMeans)

    df<-data.frame()
    df<-cbind(tmp[[1]])
    for(i in 2:length(tmp)){
      df<-cbind(df,tmp[[i]])
    }

    colnames(df)<-c(names(tmp))
    return(df)
      # OPC=as.data.frame(GetAssayData(subset(seuratobj,CellType=="OPC")
      # Pericytes=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Per")

  }
  celltype_indiv=function(celltype_means_list,
    celltype_names=c("Oligodendrocytes","Astrocytes",
    "Excitatory_Neurons","Inhibitory_Neurons","Microglia")){

    # Given a list where every element is a dataframe containing the
    # mean expression values per cell type (Gene x cell type), returns a
    # list where every element is a dataframe per cell type (gene x individual)

    commongenes<-Reduce(intersect,lapply(celltype_means_list,row.names))
    tmp<-vector("list",length=length(celltype_names))
    names(tmp)<-c(celltype_names)

    for(i in 1:length(tmp)){
      tmp[[i]]<-c(celltype_means_list[[1]][commongenes,i])
    }

    for(i in 2:length(celltype_means_list)){
      tmp[[1]]<-cbind(tmp[[1]],celltype_means_list[[i]][commongenes,1])
      tmp[[2]]<-cbind(tmp[[2]],celltype_means_list[[i]][commongenes,2])
      tmp[[3]]<-cbind(tmp[[3]],celltype_means_list[[i]][commongenes,3])
      tmp[[4]]<-cbind(tmp[[4]],celltype_means_list[[i]][commongenes,4])
      tmp[[5]]<-cbind(tmp[[5]],celltype_means_list[[i]][commongenes,5])
    }
    return(tmp)
  }

  meanslist<-lapply(splitlist,celltype_means_SCTransformed)
  pseudobulk_list<-celltype_indiv(meanslist)
  pseudobulk_list<-lapply(pseudobulk_list,function(x){colnames(x)<-samplenames
  return(x)})
  return(pseudobulk_list)

}

pseudobulk_results<-pseudobulk_scex(inputfiles)
