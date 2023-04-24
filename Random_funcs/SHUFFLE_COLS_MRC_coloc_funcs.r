



# preprocess_scex=function(inputfiles,outdir,cov.adj=TRUE,
#   SCTnormalize=TRUE,
#   CreateSeurat=FALSE,
#   covariates=c("Age","PMI","Gender")){
#
#   # inputfiles<-list.files(inputdir,pattern="dgCMatrix.rds",full.names=T)
#
#   if(CreateSeurat==TRUE){
#     scexlist<-lapply(inputfiles,readRDS)
#     seuratlist<-lapply(scexlist,CreateSeuratObject)
#
#   } else {
#     seuratlist<-lapply(inputfiles,readRDS)
#   }
#
#
#   mergeseurat<-seuratlist[[1]]
#
#   for(i in 2:length(seuratlist)){
#     mergeseurat<-merge(mergeseurat,seuratlist[[i]])
#   }
#
#  mergeseurat$CellType<-mergeseurat$ClusterIdent
#
#
#
#   # add covariates
#
#   if(cov.adj==TRUE && SCTnormalize==TRUE){
#
#     options(warn=-1)
#     message("Beginning SCTransform normalization and covariate adjustment...")
#
#
#     mergeseurat<-SCTransform(mergeseurat,variable.features.n=nrow(mergeseurat),vars.to.regress=c("Age","PMI","Gender"),verbose=FALSE)
#
#     splitlist<-SplitObject(mergeseurat,split.by="Label")
#     samplenames<-names(splitlist)
#
#
#
#     for(i in 1:length(samplenames)){
#       saveRDS(splitlist[[i]],paste0(outdir,samplenames[i],"_covadjusted_SCTransformedlabelled.rds"))
#     }
#
#
#     message("All samples SCTransform normalised, and adjusted for covariates.")
#
#
#   } else if(cov.adj==FALSE && SCTnormalize==TRUE){
#
#     options(warn=-1)
#     mergeseurat<-SCTransform(mergeseurat,variable.features.n=nrow(mergeseurat),verbose=FALSE)
#     splitlist<-SplitObject(mergeseurat,split.by="orig.ident")
#     samplenames<-names(splitlist)
#
#     for(i in 1:length(samplenames)){
#       saveRDS(splitlist[[i]],paste0(outdir,samplenames[i],"_SCTransformedlabelled.rds"))
#     }
#
#     message("SCTransform completed. No covariate adjustment done.")
#   }
#
#
# }
#
# pseudobulk_scex=function(method=c("Merged","Indiv"),inputfiles){
#
#   # This function assumes:
#
#   # 1. Input are dgcMatrix.rds files
#   # 2. Cell types are labelled
#   # 3. Covariates have been added in the meta data
#   # 4. You have done all pre-processing steps (doubletfinder,emptydrops, MT filter)
#
#   if(method=="Indiv"){
#     message("Reading in files.")
#
#     scexlist<-lapply(inputfiles,readRDS)
#     seuratlist<-lapply(scexlist,CreateSeuratObject)
#
#     mergeseurat<-seuratlist[[1]]
#
#     # SCTransform is done on the merged Seurat object
#     message("merging seurat objects..")
#     for(i in 2:length(seuratlist)){
#       mergeseurat<-merge(mergeseurat,seuratlist[[i]])
#     }
#
#
#   } else {
#     message("Reading in merged seurat object..")
#     mergeseurat<-readRDS(inputfiles)s
#   }
#
#   message("Splitting into cell types..")
#   celltypelist<-SplitObject(mergeseurat,split.by="CellType")
#   rm(mergeseurat)
#
#   options(warn=-1)
#   for(i in 1:length(celltypelist)){
#     celltypelist[[i]]<-SCTransform(celltypelist[[i]],variable.features.n=nrow(mergeseurat),vars.to.regress=c("Age","PMI","Gender"),verbose=FALSE)
#     message(paste0("SCTransformed performed for..",names(celltypelist[i])))
#
#   }
#
#   #   return(celltypelist)
#   # }
#
#   message("Obtaining expression values..")
#
#   celltypelist<-lapply(celltypelist,function(x){
#     x<-as.data.frame(GetAssayData(x,assay="SCT",slot="data"))
#   })
#
#   # options(warn=-1)
#   # message("Running SCTransform..")
#   # mergeseurat<-SCTransform(mergeseurat,variable.features.n=nrow(mergeseurat),vars.to.regress=c("Age","PMI","Gender"),verbose=FALSE)
#   # splitlist<-SplitObject(mergeseurat,split.by="orig.ident")
#   # samplenames<-names(splitlist)
#
#
#   celltype_means_SCTransformed=function(seuratobj){
#     # gets the mean per gene, per cell type, per individual (although this is done on
#     # the Pearson Residuals from the SCTransformed data)
#     library(Seurat)
#     tmp=list(oligos=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Oligo")),assay='SCT',slot='data'),
#     astros=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Astro")),assay='SCT',slot='data'),
#     Exc_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Excitatory")),assay='SCT',slot='data'),
#     Inh_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Inhibitory")),assay='SCT',slot='data'),
#     microglia=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Microglia")),assay='SCT',slot='data'),
#     OPC=as.data.frame(GetAssayData(subset(seuratobj,CellType=="OPC")),assay='SCT',slot='data'),
#     Vascular=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Vascular")),assay='SCT',slot='data'))
#
#
#     tmp=lapply(tmp,rowMeans)
#
#     df<-data.frame()
#     df<-cbind(tmp[[1]])
#     for(i in 2:length(tmp)){
#       df<-cbind(df,tmp[[i]])
#     }
#
#     colnames(df)<-c(names(tmp))
#     return(df)
#       # OPC=as.data.frame(GetAssayData(subset(seuratobj,CellType=="OPC")
#       # Pericytes=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Per")
#
#   }
#   celltype_indiv=function(celltype_means_list,
#     celltype_names=c("Oligodendrocytes","Astrocytes",
#     "Excitatory_Neurons","Inhibitory_Neurons","Microglia","OPC","Vascular")){
#
#     # Given a list where every element is a dataframe containing the
#     # mean expression values per cell type (Gene x cell type), returns a
#     # list where every element is a dataframe per cell type (gene x individual)
#
#     commongenes<-Reduce(intersect,lapply(celltype_means_list,row.names))
#     tmp<-vector("list",length=length(celltype_names))
#     names(tmp)<-c(celltype_names)
#
#     for(i in 1:length(tmp)){
#       tmp[[i]]<-c(celltype_means_list[[1]][commongenes,i])
#     }
#
#     for(i in 2:length(celltype_means_list)){
#       tmp[[1]]<-cbind(tmp[[1]],celltype_means_list[[i]][commongenes,1])
#       tmp[[2]]<-cbind(tmp[[2]],celltype_means_list[[i]][commongenes,2])
#       tmp[[3]]<-cbind(tmp[[3]],celltype_means_list[[i]][commongenes,3])
#       tmp[[4]]<-cbind(tmp[[4]],celltype_means_list[[i]][commongenes,4])
#       tmp[[5]]<-cbind(tmp[[5]],celltype_means_list[[i]][commongenes,5])
#       tmp[[6]]<-cbind(tmp[[6]],celltype_means_list[[i]][commongenes,6])
#       tmp[[7]]<-cbind(tmp[[7]],celltype_means_list[[i]][commongenes,7])
#     }
#     return(tmp)
#   }
#   message("Pseudo-bulking..")
#
#   meanslist<-lapply(splitlist,celltype_means_SCTransformed)
#   pseudobulk_list<-celltype_indiv(meanslist)
#   pseudobulk_list<-lapply(pseudobulk_list,function(x){colnames(x)<-samplenames
#   return(x)})
#   return(pseudobulk_list)
#
# }


SCT_scex=function(method=c("Merged","Indiv"),inputfiles,covs=c("Age","PMI","Gender","Sample_Source")){

  # This function assumes:

  # 1. Input are dgcMatrix.rds files
  # 2. Cell types are labelled
  # 3. Covariates have been added in the meta data
  # 4. You have done all pre-processing steps (doubletfinder,emptydrops, MT filter)

  if(method=="Indiv"){
    message("Reading in files.")

    scexlist<-lapply(inputfiles,readRDS)
    seuratlist<-lapply(scexlist,CreateSeuratObject)

    mergeseurat<-seuratlist[[1]]

    # SCTransform is done on the merged Seurat object
    message("merging seurat objects..")
    for(i in 2:length(seuratlist)){
      mergeseurat<-merge(mergeseurat,seuratlist[[i]])
    }


  } else {
    message("Reading in merged seurat object..")
    mergeseurat<-readRDS(inputfiles)
  }

  message("Splitting into cell types..")
  celltypelist<-SplitObject(mergeseurat,split.by="CellType")
  rm(mergeseurat)

  options(warn=-1)
    message(paste0("Starting SCTransform for each cell type, using following covariates: ",paste(covs,collapse=" ")))
  for(i in 1:length(celltypelist)){
    celltypelist[[i]]<-SCTransform(celltypelist[[i]],variable.features.n=nrow(celltypelist[[i]]),vars.to.regress=covs,verbose=FALSE,min_cells=5)
    message(paste0("SCTransformed performed for..",names(celltypelist[i])))

  }

  return(celltypelist)
}



pseudobulk_sct_cell=function(sct_allcells){



  tmp<-sct_allcells

  for(i in 1:length(tmp)){
    cellname=names(tmp[i])
    message(paste0("Pseudobulking ",cellname,".."))
    x<-SplitObject(tmp[[i]],split.by="Sample_ID")
    x<-lapply(x,function(x){
      x<-as.data.frame(GetAssayData(x,assay="SCT",slot="data"))
      x<-rowMeans(x)
      return(x)
    })
    x<-as.data.frame(do.call(cbind,x))
    write.table(as.matrix(x),paste0(cellname,"_pseudobulk_pseudobulk.csv"))
  }

}

labelcelltypes=function(matrix,inputdir,outdir,SCTnormalize=FALSE,LogNormalize=FALSE){
  # Given a list of known cell types/barcode, annotates Seurat object (and normalizes if specified)

  library(Seurat)
  dgcmatrix<-readRDS(paste0(inputdir,matrix))
  seuratobj<-CreateSeuratObject(counts=dgcmatrix)

  # bring in cell type labels
  labelsdir="/rds/general/user/ah3918/home/JAN21_WORK_DIR/Data/PD_Data/RAW.NONNORMALISED.COUNTS/"
  celltypelabels<-read.csv(paste0(labelsdir,"All.Samples.Cell.Labels.csv"),sep="\t")

  # needs to be fixed so as to separate into two clumns
  fixed_labels<-data.frame(do.call('rbind',strsplit(as.character(celltypelabels$x)," ",fixed=TRUE)))

  # assign cell types based on sample ID
  name<-sub(".OE.NonNormalised.filteredcounts.dgCMatrix.rds","",matrix)
  greptest<-fixed_labels[grep(name,fixed_labels$X1),]
  rownames(greptest)<-greptest$X1
  greptest<-greptest[-c(1)]
  seuratobj<-AddMetaData(object=seuratobj,metadata=greptest,col.name="CellType")

  # Normalization options
  options(warn=-1)

  if(LogNormalize==TRUE && SCTnormalize==FALSE){
    message(paste0("Performing LogNormalization and scaling for ",name,"..."))

    seuratobj<-NormalizeData(seuratobj,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)

    outfile=paste0(name,"_LogNormalizedlabelled.rds")
    saveRDS(seuratobj,paste0(outdir,outfile))
    message(paste0(name," successfully labelled, lognormalized and saved"))
  } else if(SCTnormalize==TRUE && LogNormalize==FALSE){
    message(paste0("Performing normalization by SCTransform for ",name,"..."))
    seuratobj<-SCTransform(seuratobj, verbose=FALSE)

    outfile=paste0(name,"_SCTransformedlabelled.rds")
    saveRDS(seuratobj,paste0(outdir,outfile))
    message(paste0(name," successfully labelled, SCTransform normalized and saved"))
  } else if(SCTnormalize==FALSE && LogNormalize==FALSE){
    message(paste0("Labels only for ",name,"..."))

    outfile=paste0(name,"_labelled.rds")
    saveRDS(seuratobj,paste0(outdir,outfile))
    message(paste0(name," successfully labelled and saved"))
  }
}

celltype_means=function(seuratobj){
  #gets the mean per gene, per cell type, per individual

  library(Seurat)
  tmp=list(oligos=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Oligo"))),
  astros=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Astro"))),
  Exc_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Excitatory"))),
  Inh_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Inhibitory"))),
  microglia=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Microglia"))))

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

celltype_means_SCTransformed=function(seuratobj){
  # gets the mean per gene, per cell type, per individual (although this is done on
  # the Pearson Residuals from the SCTransformed data)
  library(Seurat)
  tmp=list(oligos=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Oligo"),assay='SCT',slot='data')),
  astros=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Astro"),assay='SCT',slot='data')),
  Exc_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Excitatory"),assay='SCT',slot='data')),
  Inh_Neur=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Inhibitory"),assay='SCT',slot='data')),
  microglia=as.data.frame(GetAssayData(subset(seuratobj,CellType=="Microglia"),assay='SCT',slot='data')))

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

# --------------------------------------------------
# Get gene list from chromosome regions
# --------------------------------------------------
# inputdir="/rds/general/user/ah3918/home/JAN21_WORK_DIR/COLOC/"
# inputfile=paste0(inputdir,"genomic_coords.csv")

Genes_overlapping_coords=function(coordinatelist){
  tmp<-t(read.csv(coordinatelist,sep=",",header=F))
  tmp<-as.list(tmp[,1])

  # Create GRanges list from input ranges

  mycoords.gr = lapply(tmp, function (x) {res=strsplit(x, ':')}) %>%
     unlist %>%
     as.numeric %>%
     matrix(ncol=3, byrow=T) %>%
     as.data.frame %>%
     select(chrom=V1, start=V2, end=V3) %>%
     mutate(chrom=paste0('chr', chrom)) %>%
     makeGRangesFromDataFrame

  # Overlap with TxDb to get a list of the genes they cover
  subset<-subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), mycoords.gr)
  subset$gene_id<-as.numeric(subset$gene_id)

  # get a ist of matching genes to ID (last column of the GRanges object)

  allgenes<-as.data.frame(org.Hs.egSYMBOL)
  filteredgenelist<-allgenes[allgenes$gene_id %in% subset$gene_id,]$symbol
  return(filteredgenelist)

}


filter_df_genes=function(gene_df,genestokeep){
  # gene_matrix is a dataframe genes x sample(indiv or cell)
  # genestokeep is a character vectir of the genes to retain in df
  filteredtmp<-gene_df[row.names(gene_df) %in% genestokeep,]
  return(filteredtmp)
}


get_gene_locations=function(exp_mat){

  yourgenes<-as.character(rownames(exp_mat))

  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)

  # Create df and granges of all genes. These differ actually
  allgenes<-as.data.frame(org.Hs.egSYMBOL)
  allgranges<-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

  # get gene ids for all your genes
  gene_ids<-allgenes$gene_id[allgenes$symbol %in% yourgenes]

  # subset allgenes with the gene_ids you got above
  finalgenes<-allgenes[allgenes$symbol %in% yourgenes,]



  # the dims of granges and final genes are different.
  filteredgranges<-allgranges[allgranges$gene_id %in% gene_ids]
  df<-data.frame(filteredgranges)

  # order by gene id
  df$gene_id<-as.numeric(df$gene_id)
  df<-df[order(df$gene_id),]


  subsetgenes<-finalgenes[which(finalgenes$gene_id %in% df$gene_id),]
  df<-cbind(df,subsetgenes$gene_id,subsetgenes$symbol)

  chromlocations<-data.frame(geneid=df[,8],chr=df$seqnames,s1=df$start,s2=df$end)

  chromlocations<-chromlocations[order(chromlocations$chr),]
  colnames(chromlocations)<-c("geneid","chr","left","right")
  return(chromlocations)
  # write.table(chromlocations,"gene_locations.csv")
  #
  # message("Gene locations determined. Saved under 'gene_locations.csv' in current directory.")


}

# --------------------------------------------------
# Genotype functions
# --------------------------------------------------



# This function below will generate the genotype matrix and the locations of each snp.
# Input needs to be a list of vcfs (although if you have gone through the preprocess, just
# your outputdir will be fine. The function will automatically take the output generated from the
# preprocess step). You can also choose the minimum allele frequency, whether to take sex chromosomes
# into account, and whether you want to filter your SNPs by a list of regions, or a supplied list of SNPs.
#
# The 'Samplenamelist' is a text file from a character vector used to filter out SNPs and generate files
# only based on these samples.
#
# The snp_list is a character vector, each element is an rsid.
# The genomic_regions must be a data frame like this;
#
# chrom    start    end
# 1       120023   126000
# 2 etc .................



get_genotype_matrix=function(vcfs,
  inputfile,
  # samplenamelist,
  outdir,
  preprocess=FALSE,
  autosomalonly=TRUE,
  minmaf=0.05,
  genomewide=FALSE,
  snps_overlap=c("snp_list","genomic_regions"),
  snp_list,
  genomic_regions){

  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  library(org.Hs.eg.db)


  if(preprocess==TRUE){
    seqVCF2GDS(vcfs, paste0(outdir,"merge_test_seqArray.gds"))
  }

  genofile<<-seqOpen(inputfile)

  # 1. VERY important that the seqSetFilters are kept in the order here.
  # selecting SNPS is based on its index (its original position in the file,
  # which gets changed whenever you set a new filter. Keep this order:
  # Filter chrom -> filter snps -> filter MAF )


  if(autosomalonly==TRUE){
    seqSetFilterChrom(genofile,1:22)
  }


  annot<-seqGetData(genofile,"annotation/id")
  snp.id<-seqGetData(genofile, "variant.id")

  if(genomewide==TRUE){
    seqSetFilterCond(genofile,maf=minmaf)
    # seqSetFilter(genofile,sample.id=samplenamelist)

    annot<-seqGetData(genofile,"annotation/id")
    sample.id<-seqGetData(genofile,"sample.id")

    geno_mat<-t(seqGetData(genofile,"$dosage"))

    position<-seqGetData(genofile,"position")
    chrom<-seqGetData(genofile,"chromosome")
    chrompos_mat<-data.frame(annot,chrom,position)
    chrompos_mat<-chrompos_mat[!duplicated(chrompos_mat$annot),]

    colnames(geno_mat)<-sample.id
    geno_mat<-as.data.frame(geno_mat)
    geno_mat$rsid<-annot
    geno_mat<-geno_mat[!duplicated(geno_mat$rsid),]
    rownames(geno_mat)<-geno_mat$rsid
    geno_mat$rsid<-NULL

    chrompos_mat$chrom<-paste0("chr",chrompos_mat$chrom)
    snpnumber<-length(rownames(geno_mat))

    message(paste0("A total of ",snpnumber," snps were kept."))

    write.table(geno_mat,"genotype_012mat.csv")
    write.table(chrompos_mat,"snp_chromlocations.csv")

    seqClose(genofile)
  } else {

    if(snps_overlap=="snp_list"){

      snp_list<-readLines(snp_list)
      ids<<-which(annot %in% snp_list)

    } else if(snps_overlap=="genomic_regions"){

      position<-seqGetData(genofile,"position")
      chrom<-seqGetData(genofile,"chromosome")
      chrompos_mat<-data.frame(annot,chrom,position)

      chromGR<-makeGRangesFromDataFrame(chrompos_mat,
      seqnames.field="chrom",
      start.field="position",
      end.field="position")

      chromGR$rsid<-chrompos_mat$annot

      genomic_regions<-read.table(genomic_regions)
      overlapGR<-makeGRangesFromDataFrame(genomic_regions,
      seqnames="chrom",
      start.field="start",
      end.field="end")

      test<-subsetByOverlaps(chromGR,overlapGR)
      rsids<-as.character(as.data.frame(test)$rsid)

      ids<<-which(annot %in% rsids)
    }


    seqSetFilter(genofile,variant.id=ids)
    seqSetFilterCond(genofile,maf=minmaf)

    # seqSetFilter(genofile,sample.id=samplenamelist)

    annot<-seqGetData(genofile,"annotation/id")
    sample.id<-seqGetData(genofile,"sample.id")

    geno_mat<-t(seqGetData(genofile,"$dosage"))

    position<-seqGetData(genofile,"position")
    chrom<-seqGetData(genofile,"chromosome")
    chrompos_mat<-data.frame(annot,chrom,position)

    if(snps_overlap=="snp_list"){
      rownames(geno_mat)<-c(annot[which(annot %in% snp_list)])
    } else if(snps_overlap=="genomic_regions"){
      rownames(geno_mat)<-c(annot[which(annot %in% rsids)])
    }


    colnames(geno_mat)<-sample.id
    chrompos_mat$chrom<-paste0("chr",chrompos_mat$chrom)
    snpnumber<-length(rownames(geno_mat))

    message(paste0("A total of ",snpnumber," snps were kept."))

    write.table(geno_mat,"genotype_012mat.csv")
    write.table(chrompos_mat,"snp_chromlocations.csv")

    seqClose(genofile)
  }

}

# inputdir="/rds/general/user/ah3918/home/JAN21_WORK_DIR/PD_vcf/"
# sampleconversion<-read.csv(paste0(inputdir,"sample_ID_sample_name.csv"))
#
# PD<-sampleconversion[which(sampleconversion$Note=="Yes"),]
# PD<-PD[-which(PD$Condition=="Control"),]
# conversiontable<-data.frame(PD$sample.id,PD$sampleName)
# RES<-as.character(conversiontable$PD.sample.id)
# get_genotype_matrix(samplenamelist=RES,outdir=inputdir)

# library(GenomicRanges)
# library(GenomicFeatures)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
#
#
# gwasdir<-"/rds/general/user/ah3918/home/JAN21_WORK_DIR/GWAS_Data/"
# # PDdir<-"/rds/general/user/ah3918/home/JAN21_WORK_DIR/PD_vcf/"
# #
# pdlocs<-read.table(paste0(PDdir,"snp_chromlocations.csv"))
# #
# test<-read.table(paste0(gwasdir,"Signif_Nalls2019.csv"))
# chromlocs<-as.character(test$SNP)


chrpos_to_rsid=function(chromlocs){

  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)

  # chromlocs is a character vector, every element is
  # "chr1:10000","chr2:1201313" etc...

  chromlocs<-strsplit(chromlocs,":")
  chromlocs<-do.call(rbind.data.frame,chromlocs)
  colnames(chromlocs)<-c("Chr","Pos")

  chronly<-strsplit(as.character(chromlocs$Chr),"r")
  chronly<-do.call(rbind.data.frame,chronly)

  final<-data.frame(chrom=chronly[,2],position=chromlocs$Pos)
  final<-final[order(final$chrom),]
  final$paste<-paste0(final$chrom,":",final$position)



  grSNPS<-makeGRangesFromDataFrame(final,
  seqnames.field="chrom",
  start.field="position",
  end.field="position")



  snp<-SNPlocs.Hsapiens.dbSNP151.GRCh38
  rsids<-snpsByOverlaps(snp,grSNPS)
  rsids<-as.data.frame(rsids)$RefSNP_id

  writeLines(rsids,"Rsid_fromGWAS.txt")
}


# --------------------------------------------------
# MatrixEQTL
# # --------------------------------------------------
# exp.matrix<-"/rds/general/user/ah3918/home/JAN21_WORK_DIR/COLOC/PD_cases_only/CellTypeMeans/Oligodendrocytes_meanexpression.csv"
# genotype_matrix<-"/rds/general/user/ah3918/home/JAN21_WORK_DIR/PD_vcf/genotype_012mat.csv"
# snp_chromlocations<-"/rds/general/user/ah3918/home/JAN21_WORK_DIR/PD_vcf/snp_chromlocations.csv"
# gene_locations<-"/rds/general/user/ah3918/home/JAN21_WORK_DIR/COLOC/PD_cases_only/gene_locations.csv"
#
# calculate_ciseqtl(exp.matrix,genotype_matrix,gene_locations,snp_chromlocations,name="Oligodendrocytes")

calculate_ciseqtl=function(exp_mat,
  geno_mat,
  exp_loc,
  geno_loc,
  name,
  cisDist=1e6,
  covadj=FALSE,
  pvOutputThreshold=2e-5,
  pvOutputThreshold_cis=5e-2,
  treeqtl_threshold=0.1,
  treeqtl_only=FALSE){

  library(TreeQTL)



  exp_mat<-exp_mat
  geno_mat<-geno_mat
  exp_locs<-exp_loc
  geno_loc<-geno_loc



  # colnames(geno_mat)<-colnames(exp_mat)


  # this is in case the number of genes differ
  commongenes<-intersect(exp_locs$geneid,row.names(exp_mat))
  exp_mat<-exp_mat[commongenes,]

  # from MatrixEQTL manual: "The order of genes in gene locations
  # does not have to match the order of gene expression"

  library(MatrixEQTL)

  cvrt=SlicedData$new();
  if(covadj==TRUE){

    covariates<-read.csv("/rds/general/user/ah3918/home/JAN21_WORK_DIR/Data/PD_Data/RAW.NONNORMALISED.COUNTS/MetaData.All_Samples1.csv",sep=",",row.names=1)

    tmp<-covariates[!duplicated(covariates$IndID),]
    tmp<-tmp[which(as.character(tmp$IndID) %in% colnames(exp_mat)),]
    tmp2<-data.frame(Gender=tmp$Gender,AgeAtDeath=tmp$AgeAtDeath,PMI=tmp$PMI)
    tmp2<-t(tmp2)
    colnames(tmp2)<-as.character(tmp$IndID)
    tmp2<-as.matrix(tmp2)


    cvrt$CreateFromMatrix(as.matrix(tmp2))
  }

  snps=SlicedData$new();
  snps$CreateFromMatrix(as.matrix(geno_mat))

  gene=SlicedData$new();
  gene$CreateFromMatrix(as.matrix(exp_mat))

  if(treeqtl_only==FALSE){

    message("Performing MatrixEQTL analysis..")
    me<-Matrix_eQTL_main(snps,
                    gene,
                    cvrt = cvrt,
                    pvOutputThreshold = pvOutputThreshold,
                    useModel = modelLINEAR,
                    errorCovariance = numeric(),
                    verbose = FALSE,
                    output_file_name=paste0(name,"_trans_MatrixEQTLout.txt"),
                    output_file_name.cis =paste0(name,"_cis_MatrixEQTLout.txt"),
                    pvOutputThreshold.cis =pvOutputThreshold_cis,
                    snpspos = geno_loc,
                    genepos = exp_locs,
                    cisDist = cisDist,
                    pvalue.hist = FALSE,
                    min.pv.by.genesnp = FALSE,
                    noFDRsaveMemory = FALSE)

    message(paste0("MatrixEQTL calculated for ",name, "."))
  } else {
    message("MatrixEQTL step skipped. Loading pre-computed MatrixEQTL outputs")
  }

  message(paste0("Running TreeQTL for ",name,"..."))

  mEQTL_out_cis <- paste0(name,"_cis_MatrixEQTLout.txt")
  mEQTL_out_trans <- paste0(name,"_trans_MatrixEQTLout.txt")

  dist<-cisDist
  gene_map<-exp_locs
  gene_map$geneid<-as.character(gene_map$geneid)
  gene_map$chr<-as.character(gene_map$chr)
  gene_map$left<-as.numeric(gene_map$left)
  gene_map$right<-as.numeric(gene_map$right)



  snp_map<-geno_loc
  snp_map$annot<-as.character(snp_map$annot)
  snp_map$chrom<-as.character(snp_map$chrom)
  snp_map$position<-as.numeric(snp_map$position)

  n_tests_per_SNP <- get_n_tests_per_SNP(snp_map, gene_map, nearby = TRUE, dist = dist)
  eSNPs <- get_eSNPs(n_tests_per_SNP, mEQTL_out_cis,level1=treeqtl_threshold,level2=treeqtl_threshold)
  get_eAssociations(eSNPs, n_tests_per_SNP, mEQTL_out_cis, paste0(name,"_TreeQTL_eAssoc_cis_by_snp.txt"), by_snp = TRUE)

  # n_tests_per_gene <- get_n_tests_per_gene(snp_map, gene_map, nearby = TRUE, dist = dist)
  # eGenes <- get_eGenes(n_tests_per_gene, mEQTL_out_cis,level1=treeqtl_threshold,level2=treeqtl_threshold)
  # get_eAssociations(eGenes, n_tests_per_gene, mEQTL_out_cis, paste0(name,"_TreeQTL_eAssoc_cis_by_gene.txt"), by_snp = FALSE)


  message(paste0("Cis and Trans-eQTLs calculated for ",name))
  message(paste0("TreeQTL cutoffs determined. Saved in current dir"))

}



filter_pseudobulk=function(exp_mat,minimum_indivs){


  exp_mat<-exp_mat[rowSums(exp_mat>0)>=minimum_indivs,]
  return(exp_mat)



}


filter_genotype_matrix=function(genotypemat,filter_type=c("2.2","2.3"),nsamples=19){

  tmp<-genotypemat

  message(paste0("Filtering genotype matrix. Starting with ",nrow(genotypemat)," snps."))

  tmp$counter_0<-rowSums(tmp[1:nsamples]==0)
  tmp$counter_1<-rowSums(tmp[1:nsamples]==1)
  tmp$counter_2<-rowSums(tmp[1:nsamples]==2)

  if(filter_type=="2.2"){

    message("Retaining snps with at least 2 individuals in each of the 2 genotypic categories.")

    option1=tmp[tmp$counter_0>=2 & tmp$counter_1>=2,]
    option2=tmp[tmp$counter_0>=2 & tmp$counter_2>=2,]
    option3=tmp[tmp$counter_1>=2 & tmp$counter_2>=2,]

    intersected_rows<-unique(c(rownames(option1),rownames(option2),rownames(option3)))

  } else if (filter_type=="2.3"){

    message("Retaining snps with at least 2 individuals in each of the 3 genotypic categories.")

    tmp<-tmp[tmp$counter_0>=2 & tmp$counter_1>=2 & tmp$counter_2>=2,]

    intersected_rows<-row.names(tmp)


  }

  genotypemat<-genotypemat[intersected_rows,]

  message(paste0("Filtering complete. ",nrow(genotypemat), " snps retained."))

  return(genotypemat)
}



create_random_coords=function(numranges){

  library(TxDb.Hsapiens.UCSC.hg38.knownGene)

  rr=data.frame(chr=sample(seqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene)[1:24],numranges, replace=T))
  rr$start = apply(rr, 1, function(x) {
         round(runif(1, 0, seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)[x][[1]]), 0) }
         )
  rr$end = rr$start + runif(1, 1, 1000)  # random intervals of 1-1000 bases
  rr$strand = sample(c("+", "-"))

  rr.gr = makeGRangesFromDataFrame(rr)
  rr.gr<-as.data.frame(rr.gr)
  num<-strsplit(as.character(rr.gr$seqnames),"chr")
  num<-do.call(rbind,num)
  rr.gr$chr<-num[,2]

  df<-data.frame("chrom"=rr.gr$chr,"start"=rr.gr$start,"end"=rr.gr$end)
  return(df)




}

treeqtl_plots=function(treeqtls){
  treeqtls<-list.files(pattern="TreeQTL")
  names<-gsub("_TreeQTL_eAssoc_cis_by_snp.txt","",treeqtls)
  treeqtls<-lapply(treeqtls,read.table)

    snpgenepairlists<-lapply(treeqtls,function(x){
      x$SNP_genepair<-paste0(x$SNP,"_",x$gene)
      return(x$SNP_genepair)
  })
  names(snpgenepairlists)<-names

  library(UpSetR)
  upset(fromList(snpgenepairlists))
}

agg_cpm=function(seuratlist){
  agg_count_list<-lapply(seuratlist,function(x){
    DefaultAssay(x)<-"RNA"
    metadata<-x[[]]

    counts<-as.data.frame(GetAssayData(x))
    agg<-Matrix.utils::aggregate.Matrix(t(counts), groupings = metadata[colnames(counts),"Sample_ID"], fun = "sum")
    return(t(agg))
  })

  # now normalise
  agg_count_list<-lapply(agg_count_list,function(x){
    norm_x<-log2(cpm(x)+1)
    return(norm_x)
  })

  names(agg_count_list)<-names(seuratlist)

  #now save matrices
  for(i in 1:length(agg_count_list)){
    write.table(agg_count_list[[i]],paste0(names(agg_count_list[i]),"_agg_cpm_pseudobulk.csv"))
  }

}



agg_cpm_single=function(seuratobj){

    x<-seuratobj
    DefaultAssay(x)<-"RNA"
    metadata<-x[[]]

    counts<-as.data.frame(GetAssayData(x))
    agg<-Matrix.utils::aggregate.Matrix(t(counts), groupings = metadata[colnames(counts),"Sample_ID"], fun = "sum")
    agg<-t(agg)

    norm_x<-log2(cpm(x)+1)
  # now normalise

    return(norm_x)
}

MRC_prep_coloc=function(Step1=TRUE,
  project_name,
  inputfiles,
  outputdir,
  method=c("Merged","Indiv"),
  preprocess=TRUE,
  pseudobulk_method=c("agg_cpm","SCTransform"),
  # cov.adj_SCT=TRUE,
  # SCTnormalize=TRUE,
  #
  # covariates=c("Age","PMI","Gender"),
  Step2=TRUE,
  # samplenamelist,
  vcfs,
  gds_inputfile,
  overlap=c("genomic_regions","snp_list"),
  covs=c("Age","PMI","Gender","Sample_Source"),
  snp_list,
  genomic_regions,
  autosomalonly=TRUE,
  maf=0.05,
  genomewide=FALSE,
  filter_genotype_matrix=FALSE,
  geno_filter_type=c("2.2","2.3"),
  geno_loc_file="snp_chromlocations.csv",
  preprocess_vcf=FALSE,
  reorder_names=FALSE,
  filter_pseudobulk=TRUE,
  filter_pseudobulk_thresh=3,
  specify_genes=FALSE,
  specify_genes_list,
  # specifygenelist is a character vector
  pval_thresh_cis=5e-2,
  pval_thresh_trans=0,
  covadj_matrixeqtl=FALSE,
  treeqtl_threshold=0.1,
  cisDist=1e6,
  specify_celltypes=FALSE,
  specify_celltypes_names=NULL){

  library(GenomicRanges)
  library(GenomicFeatures)
  library(MatrixEQTL)
  library(Matrix.utils)
  library(TreeQTL)
  library(Seurat)
  library(SNPRelate)
  library(SeqArray)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  library(org.Hs.eg.db)
  library(dplyr)

  message(paste0("Starting pipeline run. Parameters used:\n covadj_matrixeqtl=",covadj_matrixeqtl,
 "\n filter_pseudobulk_thresh=",filter_pseudobulk_thresh,"\n specify_genes=",specify_genes,"\n snps used=",geno_loc_file))

  start_time<-Sys.time()
  #---------------------------------------------------
  #--------------------------------------------------
  # STEP 1:Pseudobulking
  #--------------------------------------------------
  #---------------------------------------------------
  if(Step1==TRUE){
    message("Step 1: Processing of expression matrices..")
    if(preprocess==TRUE){

      message("Pre-process=TRUE. Normalising cell types with SCTransform..")
      celltypelist<-SCT_scex(method=method,inputfiles,covs=covs)
      message("Normalisation completed, saving file.")
      saveRDS(celltypelist,paste0(outputdir,project_name,"_celltypelist_sctransformed.rds"))
    } else {
      message("Reading in pre-processed single-cell data.. ")
      celltypelist<-readRDS(paste0(outputdir,project_name,"_celltypelist_sctransformed.rds"))
    }

    message("Pseudobulking by cell type.. ")

    dir.create(paste0(outputdir,"MatrixEQTL_IO"),showWarnings=FALSE)
    matrixeqtldir<-paste0(outputdir,"MatrixEQTL_IO/")
    setwd(matrixeqtldir)

    if(pseudobulk_method=="SCTransform"){
      message("Using SCTransform normalised data and calculating mean..")
      pseudobulk_sct_cell(celltypelist)
    }else if(pseudobulk_method=="agg_cpm"){
      message("Summing counts and normalising by cpm..")
      agg_cpm(celltypelist)

    }
    message("Step 1 completed. Pseudobulk matrices written in MatrixEQTL_IO")
    message("                                    ")


  } else {
    message("Step 1 skipped. Using pre-computed pseudobulk matrices")
  }

  matrixeqtldir<-paste0(outputdir,"MatrixEQTL_IO/")
  setwd(matrixeqtldir)

  cellnames<-gsub("_pseudobulk.csv","",list.files(pattern="_pseudobulk.csv"))
  tmp<-lapply(list.files(pattern="_pseudobulk.csv"),read.table)
  names(tmp)<-cellnames

  # commongenes<-Reduce(intersect,lapply(tmp,row.names))
  # tmp<-lapply(tmp,function(x){
  #   x<-x[commongenes,]
  #   return(x)
  # })
  message("Getting gene locations for each pseudobulk matrix..")

  for(i in 1:length(tmp)){
    genelocs<-get_gene_locations(tmp[[i]])
    write.table(genelocs,paste0(names(tmp[i]),"_gene_locations.csv"))
  }

  # message("Intersection complete. Getting gene locations..")


  # get_gene_locations(tmp[[1]])






  #---------------------------------------------------
  #--------------------------------------------------
  # STEP 2: Matrix of n x p (n=sample size, p=SNPs)
  #--------------------------------------------------
  #---------------------------------------------------


  if(Step2==TRUE){

    message("Step 2: Creating genotype matrix ..")
    get_genotype_matrix(vcfs=vcfs,
      inputfile=gds_inputfile,
      preprocess=preprocess_vcf,
      outdir=outputdir,
      autosomalonly=autosomalonly,
      minmaf=maf,
      genomewide=genomewide,
      snps_overlap=overlap,
      snp_list=snp_list,
      genomic_regions=genomic_regions)
    } else {
      message("Step 2 skipped. Using pre-computed genotype matrix")
    }




  #---------------------------------------------------
  #---------------------------------------------------
  # STEP 3: Calculate MatrixEQTL
  #---------------------------------------------------
  #---------------------------------------------------

  message("Step 3: Calculating MatrixEQTL and TreeQTL..")

  setwd(matrixeqtldir)

  exp_list<-lapply(list.files(pattern="_pseudobulk.csv"),read.table)
  namelist<-gsub("_pseudobulk.csv","",list.files(pattern="_pseudobulk.csv"))
  names(exp_list)<-namelist

  # message("Matching genotype and expression matrix column names..")

  geno_mat<-read.table(list.files(pattern="genotype_012mat.csv"))
  geno_loc<-read.table(geno_loc_file)


  message("Making sure geno_loc and geno_mat match..")
  geno_mat<-geno_mat[rownames(geno_mat) %in% geno_loc$annot,]



  if(filter_genotype_matrix==TRUE){
      geno_mat<-filter_genotype_matrix(geno_mat,filter_type=geno_filter_type)
      row.names(geno_loc)<-geno_loc$annot
      geno_loc<-geno_loc[rownames(geno_mat),]
      row.names(geno_loc)<-rep(1:nrow(geno_loc))
  }

  if(filter_pseudobulk==TRUE){

    message(paste0("filter_pseudobulk=TRUE. Keeping genes with minimum ",filter_pseudobulk_thresh, " individuals."))
    exp_list<-lapply(exp_list,filter_pseudobulk,minimum_indivs=filter_pseudobulk_thresh)
  }

  if(specify_genes==TRUE){

    exp_list<-lapply(exp_list,function(x){
      return(x[rownames(x) %in% specify_genes_list,])
    })
  }

  # each pseudobulk matrix has its own gene locations file
  exp_loc<-lapply(list.files(pattern="gene_locations.csv"),read.table)

  if(specify_celltypes==TRUE){
    message("specify_celltypes=TRUE. Specific cell types selected.")
    message(paste0("Performing analysis for ",specify_celltypes_names))

    namelist<-specify_celltypes_names
    exp_list<-exp_list[namelist]
  }
  
  for(i in 1:length(exp_list)){
    message(paste0("Ordering column names for ",namelist[i]))
    geno_mat<-geno_mat[colnames(exp_list[[i]])]
    message("Shuffling colnames..")
    exp_list[[i]]<-exp_list[[i]][sample(ncol(exp_list[[i]]))]

    calculate_ciseqtl(exp_mat=exp_list[[i]],
    exp_loc=exp_loc[[i]],
    geno_mat=geno_mat,
    geno_loc=geno_loc,
    name=namelist[i],
    pvOutputThreshold= pval_thresh_trans,
    pvOutputThreshold_cis=pval_thresh_cis,
    treeqtl_threshold=treeqtl_threshold,
    covadj=covadj_matrixeqtl)
  }



  end_time<-Sys.time()

  duration=end_time - start_time

  message("Analysis complete!")
  message(paste0("Runtime: ",duration))

}
