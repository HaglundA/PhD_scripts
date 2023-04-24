




preprocess_scex=function(inputfiles,outdir,cov.adj=TRUE,SCTnormalize=TRUE,LogNormalize=FALSE){

  # inputfiles<-list.files(inputdir,pattern="dgCMatrix.rds",full.names=T)
  scexlist<-lapply(inputfiles,readRDS)
  seuratlist<-lapply(scexlist,CreateSeuratObject)

  message("Merging all samples...")
  mergeseurat<-Reduce(merge,seuratlist)



  # 1. Label cell types
  labelsdir="/rds/general/user/ah3918/home/JAN21_WORK_DIR/Data/PD_Data/RAW.NONNORMALISED.COUNTS/"
  celltypelabels<-read.csv("/mnt/data/ALEX/COLOC_FILES/All.Samples.Cell.Labels.csv",sep="\t")


  celltypelabels<-data.frame(do.call('rbind',strsplit(as.character(celltypelabels$x)," ",fixed=TRUE)))
  rownames(celltypelabels)<-celltypelabels$X1

  celltypelabels<-as.data.frame(celltypelabels[colnames(mergeseurat),])
  celltypelabels$X1<-NULL
  mergeseurat<-AddMetaData(mergeseurat,metadata=celltypelabels,col.name="CellType")



  # add covariates

  if(cov.adj==TRUE && SCTnormalize==TRUE){

    options(warn=-1)
    message("Beginning SCTransform normalization and covariate adjustment...")

    covariates<-read.csv("/mnt/data/ALEX/COLOC_FILES/MetaData.All_Samples1.csv",sep=",",row.names=1)

    covariates<-covariates[colnames(mergeseurat),]
    covariates<-covariates[order(match(colnames(mergeseurat),rownames(covariates))),c("AgeAtDeath","PMI","Gender")]

    mergeseurat[["Age"]]=covariates["AgeAtDeath"]
    mergeseurat[["PMI"]]=covariates["PMI"]
    mergeseurat[["Gender"]]=covariates["Gender"]

    # split by cell type and run SCtransform for each cell type
    message("Splitting by cell type..")
    celltypelist<-SplitObject(mergeseurat,split.by="CellType")
    for(i in 1:length(celltypelist)){
      celltypelist[[i]]<-SCTransform(celltypelist[[i]],variable.features.n=nrow(celltypelist[[i]]),vars.to.regress=c("Age","PMI","Gender"),verbose=FALSE)
    }

    message("Merging them back together...")

    mergeseurat<-Reduce(merge,celltypelist)
    splitlist<-SplitObject(mergeseurat,split.by="orig.ident")
    samplenames<-names(splitlist)



    for(i in 1:length(samplenames)){
      saveRDS(splitlist[[i]],paste0(outdir,samplenames[i],"_covadjusted_SCTransformedlabelled.rds"))
    }


    message("All samples SCTransform normalised, and adjusted for covariates.")


  } else if(cov.adj==FALSE && SCTnormalize==TRUE){

    options(warn=-1)
    mergeseurat<-SCTransform(mergeseurat,variable.features.n=nrow(mergeseurat),verbose=FALSE)
    splitlist<-SplitObject(mergeseurat,split.by="orig.ident")
    samplenames<-names(splitlist)

    for(i in 1:length(samplenames)){
      saveRDS(splitlist[[i]],paste0(outdir,samplenames[i],"_SCTransformedlabelled.rds"))
    }

    message("SCTransform completed. No covariate adjustment done.")

  } else if(cov.adj==FALSE && LogNormalize==TRUE){

    NormalizeData(mergeseurat,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
    splitlist<-SplitObject(mergeseurat,split.by="orig.ident")
    samplenames<-names(splitlist)

    for(i in 1:length(samplenames)){
      saveRDS(splitlist[[i]],paste0(outdir,samplenames[i],"_SCTransformedlabelled.rds"))
    }

    message("LogNormalization completed. No covariate adjustment done.")

  }



}
