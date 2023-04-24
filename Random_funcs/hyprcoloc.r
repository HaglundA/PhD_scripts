








# NEW COLOC PIPELINE:


select_regions=function(gwas,window,pval){


  library(parallel)
  message("Splitting into chromosomes..")
  gwas_chr_list<-split(gwas,gwas$chr)

  kb_window=window/1e3
  message("Selecting lead index snps per region..")
  gwas_clump<-mclapply(gwas_chr_list,ld_clump,clump_kb=kb_window,clump_r2=0.001,clump_p=pval,mc.cores=4)

  gwas_clump<-gwas_clump[sapply(gwas_clump, nrow)>0]
  gwas_chr_list<-gwas_chr_list[names(gwas_clump)]

  populate_clumps=function(clumped_gwas_chr,gwas_chr,window=window,pval=pval){

    signif_snps<-clumped_gwas_chr
    x<-gwas_chr

    signif_snps$pos<-as.numeric(signif_snps$pos)
    signif_snps$start<-signif_snps$pos-window
    signif_snps$end<-signif_snps$pos+window

    start<-signif_snps[1,]$start
    end<-signif_snps[1,]$end
    df<-subset(x,pos>start & pos<end)
    df$n_snps_region<-rep(nrow(df),nrow(df))
    df$signif_snp_region<-rep(signif_snps[1,]$rsid_hg38,nrow(df))


    message("Populating regions..")
    if(nrow(signif_snps)>1){
      for(i in 2:nrow(signif_snps)){
        start<-signif_snps[i,]$start
        end<-signif_snps[i,]$end
        df_2<-subset(x,pos>start & pos<end)
        df_2$n_snps_region<-rep(nrow(df_2),nrow(df_2))
        df_2$signif_snp_region<-rep(signif_snps[i,]$rsid_hg38,nrow(df_2))

        df<-rbind(df,df_2)
      }
      return(df)
    } else {
      return(df)
    }
  }


  window=window
  for(i in 1:length(gwas_chr_list)){
    test[[i]]<-populate_clumps(clumped_gwas=gwas_clump[[i]],gwas_chr=gwas_chr_list[[i]],window=window)
    }

  test<-do.call(rbind,test)
  return(test)


}

createdf=function(x,uniquegene){

  generows<-which(x$gene==as.character(uniquegene))

  df<-data.frame(SNP=as.character(x$SNP[generows]),
  beta=x$beta[generows],
  se=x$std_error[generows],
  GWAS_se=x$GWAS_se[generows],
  GWAS_beta=x$GWAS_beta[generows])
  return(df)

}

preprocess_mateqtlouts=function(meqtloutdir){
  message("Reading in MatrixEQTL outputs..")
  mateqtlouts<-list.files(meqtloutdir,pattern="MatrixEQTL",full.names=T)
  mateqtlouts<-lapply(mateqtlouts,read.table,head=T)

  cellnames<-gsub("_cis_MatrixEQTLout.txt","",list.files(meqtloutdir,pattern="MatrixEQTL"))

  names(mateqtlouts)<-cellnames


  rowsbefore<-nrow(mateqtlouts[[1]])
  message(paste0("Beginning with pre-processing. Assuming all cell types share same number of cis-EQTLs, starting with ",rowsbefore, " cis-EQTLs"))

  mateqtlouts<-lapply(mateqtlouts,function(x){
    row.names(x)<-paste0(x$SNP,"_",x$gene)
    return(x)
  })

  # hyprcoloc can't handle zero values. These are filtered out
  mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[!x$t.stat==0,]
    x<-x[!x$t.stat==-Inf,]
    x$std_error=x$beta/x$t.stat
    return(x)})

  message("Standard errors determined. Removing zero values and intersecting row names.")
  #
  # truevec<-list()
  # for(i in 1:length(mateqtlouts)){
  #     truevec[[i]]<-all.equal(mateqtlouts[[i]]$SNP,mateqtlouts[[1]]$SNP)
  # }
  #
  # if(all(truevec==TRUE)==FALSE) stop("SNP IDs don't match. Are you sure you're inputting genome-wide data? All MatrixEQTL outputs must have the same associations.")

  # If we want to run hyprcoloc between all cell types, they must all have the same cis-eqtls.
  commonrows<-Reduce(intersect,lapply(mateqtlouts,row.names))

  mateqtlouts<<-lapply(mateqtlouts,function(x){x<-x[commonrows,]
    return(x)})
  rowsdiff=rowsbefore-nrow(mateqtlouts[[1]])

  message(paste0("Intersection done. ",rowsdiff, " cis-EQTLs removed from analysis."))
  return(mateqtlouts)
}

appendGWAS=function(x,GWAS){

  x$GWAS_se<-GWAS$se[match(x$SNP,GWAS$rsid_hg38)]
  x$GWAS_beta<-GWAS$b[match(x$SNP,GWAS$rsid_hg38)]
  return(x)
}

get_genes_per_region=function(processed_gwas_region,gene_locs){

  gene_locs<-split(gene_locs,gene_locs$chr)
  gene_locs<-gene_locs[which(names(gene_locs)==unique(processed_gwas_region$chr))]
  gene_locs<-do.call(rbind,gene_locs)
  rownames(gene_locs)<-rep(1:nrow(gene_locs))

  startpos<-min(processed_gwas_region$pos)
  endpos<-max(processed_gwas_region$pos)

  gene_locs<-subset(gene_locs,left>startpos & left<endpos)

  return(gene_locs$geneid)

}

prep_coloc_data=function(mateqtlouts,processed_gwas,gene_locs){

  if(class(processed_gwas)=="data.frame"){
    processed_gwas<-split(processed_gwas,processed_gwas$signif_snp_region)
  }
  biglist<-list()
  message("Intersecting GWAS region with eQTL genes and snps..")
  for(i in 1:length(processed_gwas)){

    # message(paste0("Selecting eQTLs for region ",i))
    tmp_gwas<-processed_gwas[[i]]
    genestokeep<-get_genes_per_region(tmp_gwas,gene_locs)

    biglist[[i]]<-lapply(mateqtlouts,function(x){
      x<-x[x$gene %in% genestokeep,]
      x<-x[x$SNP %in% tmp_gwas$rsid,]
      return(x)
    })
  }
  names(biglist)<-names(processed_gwas)

  message("Generating beta and std error data frames..")
  biglist_se<-list()
  for(i in 1:length(biglist)){
    biglist_se[[i]]<-lapply(biglist[[i]],function(x){
          x<-x[c("SNP","gene","std_error")]
          x<-x%>% pivot_wider(names_from = gene, values_from = std_error)
          x<-as.data.frame(x)
          rownames(x)<-x$SNP
          x$SNP<-NULL
          return(x)})
        }

  biglist_beta<-list()
  for(i in 1:length(biglist)){
    biglist_beta[[i]]<-lapply(biglist[[i]],function(x){
      x<-x[c("SNP","gene","beta")]
      x<-x%>% pivot_wider(names_from = gene, values_from = beta)
      x<-as.data.frame(x)
      rownames(x)<-x$SNP
      x$SNP<-NULL
      return(x)})
    }


    for(i in 1:length(biglist_se)){
      biglist_se[[i]]<-do.call(cbind,biglist_se[[i]])
      biglist_beta[[i]]<-do.call(cbind,biglist_beta[[i]])
    }

    names(biglist_se)<-names(biglist)
    names(biglist_beta)<-names(biglist)
    biglist_se<-biglist_se[sapply(biglist_se, nrow)>0]
    biglist_beta<-biglist_beta[sapply(biglist_beta, nrow)>0]

    biglist_se<-lapply(biglist_se,appendGWAS_se,GWAS=processed_gwas)
    biglist_beta<-lapply(biglist_beta,appendGWAS_beta,GWAS=processed_gwas)
    biglist_beta<-lapply(biglist_beta,as.matrix)
    biglist_se<-lapply(biglist_se,as.matrix)

    return_list<-list(biglist_se,biglist_beta)
    names(return_list)<-c("Betas","StdErrors")
    return(return_list)
}

run_hyprcoloc=function(biglist_beta,biglist_se,processed_gwas){


  hyprcoloc_results<-list()
  for(i in 1:length(biglist_beta)){
        traits_list[[i]]<-colnames(biglist_beta[[i]])
        snp_id[[i]]<-rownames(biglist_beta[[i]])
        hyprcoloc_results[[i]]<-hyprcoloc(biglist_beta[[i]],
          biglist_se[[i]], trait.names=traits_list[[i]], snp.id=snp_id[[i]])
      }

  final<-lapply(hyprcoloc_results,function(x){
    x<-as.data.frame(x$results)
    x<-x[order(x$posterior_prob,decreasing=T),]
    x<-x[which(x$posterior_prob>0),]
    return(x)})

  final<-final[sapply(final, nrow)>0]
  hyprcoloc_res_2<-do.call(rbind,final)

  final<-append_gwas_info(final,processed_gwas)
  final<-final[c("traits","regional_prob","posterior_prob","posterior_explained_by_snp",
               "candidate_snp","GWAS_pval","n_snps_region","index_snp")]


  final<-final[order(final$GWAS_pval),]
  return(final)

}

# This is for the region selection - 08MAR21
# -====================================




# # create GR object from regions of interest
# overlapGR<-makeGRangesFromDataFrame(genomic_regions,
# seqnames="chrom",
# start.field="start",
# end.field="end")
#
# # create GR from snps
#
# snps<-SNPlocs.Hsapiens.dbSNP151.GRCh38
# my_snps<-as.character(mateqtlouts[[1]]$SNP)
#
# locs<-snpsById(snps,my_snps)
# locs<-as.data.frame(locs)
# locsGR<-makeGRangesFromDataFrame(locs,
# seqnames="seqnames",
# start.field="pos",
# end.field="pos")
# values(locsGR)<-DataFrame(rsid=locs$RefSNP_id)
#
# # overlap
#
# overlap<-as.data.frame(subsetByOverlaps(locsGR,overlapGR))
#
# # intersect with overlapGR
#
# mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% overlap$rsid),]})
#
# # subset matrices by the remaining snps

# -====================================
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(org.Hs.eg.db)
library(dplyr)
library(hyprcoloc)
library(stringr)

run_hyprcoloc=function(GWAS,
  meqtloutdir,
  snplocs,
  running_mode=c("Genome_wide","Regions_of_interest","snplist"),
  genomic_regions,
  snplist,
  GWASsignif=10e-6,
  preprocess=FALSE,
  plots=TRUE,
  APITOKEN="f3d054e6c0ee",
  outdir="/.",
  celltype_names=c("Astrocytes","Excitatory_Neurons","Inhibitory_Neurons","Microglia","Oligodendrocytes")){


  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  library(org.Hs.eg.db)
  library(dplyr)
  library(hyprcoloc)
  library(stringr)
  library(LDlinkR)





  # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 1. Reading in data and filtering -=-=-=-=-=-=-===
  # -============================================================================================================

  if(preprocess==TRUE){

    message("Reading in MatrixEQTL outputs..")
    mateqtlouts<-list.files(meqtloutdir,pattern="MatrixEQTL",full.names=T)
    mateqtlouts<-lapply(mateqtlouts,read.table,head=T)

    cellnames<-gsub("_cis_MatrixEQTLout.txt","",list.files(meqtloutdir,pattern="MatrixEQTL"))

    names(mateqtlouts)<-cellnames


    rowsbefore<-nrow(mateqtlouts[[1]])
    message(paste0("Beginning with pre-processing. Assuming all cell types share same number of cis-EQTLs, starting with ",rowsbefore, " cis-EQTLs"))

    mateqtlouts<-lapply(mateqtlouts,function(x){
      row.names(x)<-paste0(x$SNP,"_",x$gene)
      return(x)
    })

    # hyprcoloc can't handle zero values. These are filtered out
    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[!x$t.stat==0,]
      x<-x[!x$t.stat==-Inf,]
      x$std_error=x$beta/x$t.stat
      return(x)})

    message("Standard errors determined. Removing zero values and intersecting row names.")
    #
    # truevec<-list()
    # for(i in 1:length(mateqtlouts)){
    #     truevec[[i]]<-all.equal(mateqtlouts[[i]]$SNP,mateqtlouts[[1]]$SNP)
    # }
    #
    # if(all(truevec==TRUE)==FALSE) stop("SNP IDs don't match. Are you sure you're inputting genome-wide data? All MatrixEQTL outputs must have the same associations.")

    # If we want to run hyprcoloc between all cell types, they must all have the same cis-eqtls.
    commonrows<-Reduce(intersect,lapply(mateqtlouts,row.names))

    mateqtlouts<<-lapply(mateqtlouts,function(x){x<-x[commonrows,]
      return(x)})
    rowsdiff=rowsbefore-nrow(mateqtlouts[[1]])

    message(paste0("Intersection done. ",rowsdiff, " cis-EQTLs removed from analysis."))



    if(running_mode=="Regions_of_interest"){

      message("Selecting cis-EQTLs based on regions of interest..")

      genomic_regions<-read.table(genomic_regions)
      overlapGR<-makeGRangesFromDataFrame(genomic_regions,
      seqnames="chrom",
      start.field="start",
      end.field="end")

      snps<-SNPlocs.Hsapiens.dbSNP151.GRCh38

      # if the matrixeqtloutput is genomewide (as it should be), they should all have exactly
      # the same number of SNPs and in the same order. So it doesn't matter that I
      # just take the first element of that list.

      my_snps<-as.character(mateqtlouts[[1]]$SNP)
      snps<-SNPlocs.Hsapiens.dbSNP151.GRCh38

      locs<-snpsById(snps,my_snps)
      locs<-as.data.frame(locs)
      locsGR<-makeGRangesFromDataFrame(locs,
      seqnames="seqnames",
      start.field="pos",
      end.field="pos")
      values(locsGR)<-DataFrame(rsid=locs$RefSNP_id)

      overlap<-as.data.frame(subsetByOverlaps(locsGR,overlapGR))

      mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% overlap$rsid),]
      return(x)})
      cisremain<-nrow(mateqtlouts[[1]])

      message(paste0("Selection complete. ",cisremain," cis-EQTLs remain in the analysis"))

    } else if (running_mode=="snplist"){
      message("Intersecting with SNPs of interest..")
      snplist<-read.table(snplist)
      mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% snplist),]
      return(x)})

      cisremain<-nrow(mateqtlouts[[1]])
      message(paste0("Selection complete. ",cisremain," cis-EQTLs remain in the analysis"))
    }





  # read in GWAS data
    message("Reading in and overlapping with GWAS data..")
    GWAS<-read.table(GWAS)
    GWAS<-GWAS[GWAS$p<GWASsignif,]


    # filter meqtl matrices based on overlap with GWAS
    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% GWAS$rsid_hg38),]
      return(x)})

    # second step of intersection
    commonrows<-Reduce(intersect,lapply(mateqtlouts,row.names))
    mateqtlouts<-lapply(mateqtlouts,function(x){
      x<-x[commonrows,]
      return(x)
    })

    message("Saving pre-processed list under 'mateqtlouts.rds' in current dir. You can re-start the analysis from here by specifying preprocess=FALSE.")

    saveRDS(mateqtlouts,"mateqtlouts.rds")
  }

  if(preprocess==FALSE){
    message("Reading in pre-processed list and GWAS data..")
    mateqtlouts<-readRDS("mateqtlouts.rds")

    GWAS<-read.table(GWAS)
    GWAS<-GWAS[GWAS$p<GWASsignif,]
    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% GWAS$rsid_hg38),]
      return(x)})

  }



  # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 2. Generating the gene-specific SE and Beta matrices -=-=-=-=-=-=-===
  # -============================================================================================================


  createdf=function(x,uniquegene){

    generows<-which(x$gene==as.character(uniquegene))

    df<-data.frame(SNP=as.character(x$SNP[generows]),
    beta=x$beta[generows],
    se=x$std_error[generows],
    GWAS_se=x$GWAS_se[generows],
    GWAS_beta=x$GWAS_beta[generows])
    return(df)

  }

  appendGWAS=function(x,GWAS){

    x$GWAS_se<-GWAS$se[match(x$SNP,GWAS$rsid_hg38)]
    x$GWAS_beta<-GWAS$b[match(x$SNP,GWAS$rsid_hg38)]
    return(x)
  }

  # append GWAS_se and beta values
  mateqtlouts<-lapply(mateqtlouts,appendGWAS,GWAS=GWAS)

  uniquegenes<-unique(mateqtlouts[[1]]$gene)
  uniquegenes<-as.character(uniquegenes)

  # create a big list which will contain multiple lists (one for each cell type)
  biglist<-vector("list",length=length(mateqtlouts))

  for(i in 1:length(mateqtlouts)){
    biglist[[i]]<-lapply(uniquegenes,createdf,x=mateqtlouts[[i]])
  }

  # gene_list_astro<-lapply(uniquegenes,createdf,x=mateqtlouts[[1]],GWAS=GWAS)

  HC_input<-biglist[[1]]
  for(i in 2:length(biglist)){
    HC_input<-mapply(cbind,HC_input,biglist[[i]],SIMPLIFY=FALSE)
  }

  names(HC_input)<-uniquegenes



  # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 3. Run hyprcoloc -=-=-=-=-=-=-===
  # -============================================================================================================

  calculate_hypr=function(x,names){

    row.names(x)<-x$SNP

    # betas
    x2<-x[,which(colnames(x)=="beta")]
    colnames(x2)<-paste0(names,colnames(x2))
    x2$GWAS_beta<-x$GWAS_beta
    row.names(x2)<-x$SNP

    x2<-x2 %>% mutate_if(is.null, as.numeric)
    x2<-x2 %>% mutate_if(is.character, as.numeric)
    x3<-as.matrix(x2)

    # standard errorrs

    y2<-x[,which(colnames(x)=="se")]
    colnames(y2)<-paste0(names,"_",colnames(y2))
    y2$GWAS_se<-x$GWAS_se
    row.names(y2)<-x$SNP
    y2 <- y2 %>% mutate_if(is.null, as.numeric)
    y2 <- y2 %>% mutate_if(is.character, as.numeric)
    y3 <- as.matrix(y2)

    traits <- colnames(x3)
    variants <- rownames(x3)
    binary.traits <- as.numeric(!grepl("GWAS", traits))

    # calculate LD mat
    # if(nrow(x3)>=2 && nrow(x3)<=1000){
    #      message("SNPs within range. Calculating LD.")
    #      suppressMessages(cor<-LDmatrix(row.names(x3), pop = "CEU", r2d = "r2", token ="f3d054e6c0ee", file = FALSE))
    #      rownames(cor)<-cor$RS_number
    #      cor$RS_number<-NULL
    #      ldmat<-as.matrix(cor)
    #
    #  } else if(nrow(x3)>1000){
    #     message("More than 1000 snps, cutting matrix")
    #      x3<-x3[1:1000,]
    #      y3<-y3[1:1000,]
    #      suppressMessages(cor<-LDmatrix(row.names(x3), pop = "CEU", r2d = "r2", token ="f3d054e6c0ee", file = FALSE))
    #      rownames(cor)<-cor$RS_number
    #      cor$RS_number<-NULL
    #      ldmat<-as.matrix(cor)
    #    }



    if (nrow(x3) > 1) {
      res <- hyprcoloc(x3, y3, trait.names = traits, snp.id = variants,snpscores = TRUE)
      return(res$results)
    }

  }

  get_SNPScores=function(x,names){
   library(dplyr)
   library(hyprcoloc)
   row.names(x)<-x$SNP

   # betas
   x2<-x[,which(colnames(x)=="beta")]
   colnames(x2)<-paste0(names,colnames(x2))
   x2$GWAS_beta<-x$GWAS_beta
   row.names(x2)<-x$SNP

   x2<-x2 %>% mutate_if(is.null, as.numeric)
   x2<-x2 %>% mutate_if(is.character, as.numeric)
   x3<-as.matrix(x2)

   # standard errorrs

   y2<-x[,which(colnames(x)=="se")]
   colnames(y2)<-paste0(names,"_",colnames(y2))
   y2$GWAS_se<-x$GWAS_se
   row.names(y2)<-x$SNP
   y2 <- y2 %>% mutate_if(is.null, as.numeric)
   y2 <- y2 %>% mutate_if(is.character, as.numeric)
   y3 <- as.matrix(y2)

   traits <- colnames(x3)
   variants <- rownames(x3)
   binary.traits <- as.numeric(!grepl("GWAS", traits))

   # if(nrow(x3)>=2 && nrow(x3)<=1000){
   #        message("SNPs within range. Calculating LD.")
   #        suppressMessages(cor<-LDmatrix(row.names(x3), pop = "CEU", r2d = "r2", token ="f3d054e6c0ee", file = FALSE))
   #        rownames(cor)<-cor$RS_number
   #        cor$RS_number<-NULL
   #        ldmat<-as.matrix(cor)
   #
   #    } else if(nrow(x3)>1000){
   #       message("More than 1000 snps, cutting matrix")
   #        x3<-x3[1:1000,]
   #        y3<-y3[1:1000,]
   #        suppressMessages(cor<-LDmatrix(row.names(x3), pop = "CEU", r2d = "r2", token ="f3d054e6c0ee", file = FALSE))
   #        rownames(cor)<-cor$RS_number
   #        cor$RS_number<-NULL
   #        ldmat<-as.matrix(cor)
   #      }
   #
   if (nrow(x3) > 1) {
     res <- hyprcoloc(x3, y3, trait.names = traits, snp.id = variants, snpscores = TRUE)
     return(res$snpscores)
    }
  }

  celltype_names<-names(mateqtlouts)

  message("Running hyprcoloc..")
  results<-lapply(HC_input,calculate_hypr,names=celltype_names)
  snpscores<-lapply(HC_input,get_SNPScores,names=celltype_names)
  names(results)<-as.character(uniquegenes)
  names(snpscores)<-as.character(uniquegenes)

  finaldf<-do.call(rbind,results)
  finaldf$gene<-rownames(finaldf)
  finaldf$gene<-substr(finaldf$gene,1,nchar(finaldf$gene)-2)
  rownames(finaldf)<-c(1:nrow(finaldf))
  saveRDS(finaldf,"All_hyprcoloc_results.rds")
  saveRDS(snpscores,"SNPScores.rds")
  saveRDS(results,"results.rds")


  # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 4. Generate plots  -=-=-=-=-=-=-===
  # -============================================================================================================

  getpval=function(x,mateqtlouts){

    snpgenepair<-paste0(x$candidate_snp,"_",x$gene)

    traits<-x$traits
    traits<-strsplit(traits,", ")
    celltypename<-traits[[1]][1]
    if(length(traits[[1]])>2){
        celltypename1<-traits[[1]][1]
        celltypename2<-traits[[1]][2]

        df1<-mateqtlouts[[which(names(mateqtlouts)==celltypename1)]]
        pval1<-df1[snpgenepair,]$p.value
        df2<-mateqtlouts[[which(names(mateqtlouts)==celltypename2)]]
        pval2<-df2[snpgenepair,]$p.value

        pval<-paste0(celltypename1,": ",pval1,"\n",celltypename2,": ",pval2)


    } else {
      df<-mateqtlouts[[which(names(mateqtlouts)==celltypename)]]
      pval<-df[snpgenepair,]$p.value
      }

    return(pval)
  }

  get_BBFDR=function(x,mateqtlouts){

    snpgenepair<-paste0(x$candidate_snp,"_",x$gene)

    traits<-x$traits
    traits<-strsplit(traits,", ")
    celltypename<-traits[[1]][1]
    if(length(traits[[1]])>2){
        celltypename1<-traits[[1]][1]
        celltypename2<-traits[[1]][2]
        df1<-mateqtlouts[[which(names(mateqtlouts)==celltypename1)]]
        pval1<-df1[snpgenepair,]$BBFDR
        df2<-mateqtlouts[[which(names(mateqtlouts)==celltypename2)]]
        pval2<-df2[snpgenepair,]$BBFDR

        pval<-paste0(celltypename1,": ",pval1,"\n",celltypename2,": ",pval2)


      } else {
        df<-mateqtlouts[[which(names(mateqtlouts)==celltypename)]]
        pval<-df[snpgenepair,]$BBFDR
        }

      return(pval)
  }


  # filter based on which colocalisation have a posterior prob
  finaldf$traits<-str_remove_all(finaldf$traits,paste(c("beta","beta.1","_beta",".1",".2",".3",".4",".5",".6",".7"),collapse="|"))
  finaldf_hits<-finaldf[which(finaldf$posterior_prob>0),]

  # add in the GWAS p-value
  finaldf_hits$GWAS_pval<-cbind(GWAS[match(finaldf_hits$candidate_snp,GWAS$rsid_hg38),]$p)
  # Added in 21JUL21
  finaldf_hits$GWAS_pval<-as.numeric(finaldf_hits$GWAS_pval)


  # now the cis-eqtl pvalue
  pvals<-vector()
  for(i in 1:nrow(finaldf_hits)){
    pvals<-c(pvals,getpval(finaldf_hits[i,],mateqtlouts))
  }

  finaldf_hits$CisEQTL_pval<-cbind(pvals)

  # finally, check if available in the treeqtl output
  meqtloutdir=meqtloutdir
  treeqtlfiles=list.files(meqtloutdir,pattern="_TreeQTL_eAssoc_cis_by_snp.txt",full.names=T)
  treeqtl_names<-gsub(meqtloutdir,"",treeqtlfiles)
  treeqtl_names<-gsub("_TreeQTL_eAssoc_cis_by_snp.txt","",treeqtl_names)
  treeqtl_names<-gsub("/","",treeqtl_names)


  treeqtlfiles<-lapply(treeqtlfiles,read.table,head=T)
  names(treeqtlfiles)<-treeqtl_names

  treeqtlfiles<-lapply(treeqtlfiles,function(x){
    row.names(x)<-paste0(x$SNP,"_",x$gene)
    return(x)
  })

  pvals<-vector()
  for(i in 1:nrow(finaldf_hits)){
    pvals<-c(pvals,get_BBFDR(finaldf_hits[i,],treeqtlfiles))
  }
  finaldf_hits$TreeQTL_BBFDR<-cbind(pvals)



  write.table(finaldf,"finaldf.txt")
  write.table(finaldf_hits,"finaldf_hits.txt")
  message("Analysis complete. Co-localisation hits saved under 'finaldf_hits.txt'.")
  message("To create co-localisation plots, use the create_gassocplo2 function.")
  message("You must provide the finaldf_hits file, the snp_locations and the SNPScores file.")

}




for(i in 1:nrow(finaldf_hits)){

}





create_gassocplot2=function(finaldf_row,snplocs,genelocs,snpscores,APITOKEN="yourtoken"){
  library(LDlinkR)
  library(gassocplot2)
  if(APITOKEN=="yourtoken"){
    message("Error: Please input an API token with APITOKEN=yourtoken. This function calculates LD scores between the SNPs of interest.")
    message("It utilizes LDLinkR, which requires an API token from https://ldlink.nci.nih.gov/?tab=apiaccess")
  } else {

    snpscores<-snpscores
    locs<-snplocs
    genelocs<-genelocs
    gene=finaldf_row$gene

    library(stringr)
    scores<-snpscores[which(names(snpscores)==finaldf_row$gene)][[1]]
    if(scores[1]=="NULL"){
      scores<-snpscores[which(names(snpscores)==finaldf_row$gene)][[1]][2]
      scores<-as.data.frame(scores)
    } else {
      scores<-as.data.frame(snpscores[which(names(snpscores)==finaldf_row$gene)][[1]])
    }


    colnames(scores)<-finaldf_row$traits
    traits<-colnames(scores)

    l<-locs[locs$annot %in% row.names(scores),]
    l$chrom<-str_remove(l$chrom,"chr")
    colnames(l)<-c("marker","chr","pos")

    rownames(l)<-l$marker


    # LDMAT
    # LDLink only takes a maximum of 1000 snps.
    snps_ld<-rownames(scores)

    if(nrow(scores)>1000){
      snps_ld<-snps_ld[1:1000]
    }

    cor<-LDmatrix(snps_ld, pop = "CEU", r2d = "r2", token = APITOKEN, file = FALSE)
    rownames(cor)<-cor$RS_number
    cor$RS_number<-NULL
    scores<-as.data.frame(scores[rownames(cor),])
    rownames(scores)<-rownames(cor)
    l<-l[rownames(scores),]

    if(nrow(scores)<10){
      top3<-scores[order(scores[,1],decreasing=T),,drop=FALSE]
      top3<-top3[2:3,,drop=FALSE]
      top3<-row.names(top3)
    } else {
      top3<-NULL
    }


    # get gene locs
    leadsnp<-rownames(scores)[which(scores==max(scores))][1]
    leadsnp_x<-l[leadsnp,]$pos
    genepos_df=genelocs[which(genelocs$geneid==gene),]
    if(leadsnp_x<genepos_df$right){
      message("gene is downstream")
      xmax=genepos_df$right+200
      xmin=leadsnp_x-200
    }else{
      message("gene is upstream")
      xmin=genepos_df$left-200
      xmax=leadsnp_x+200
    }
      # Create plot
    plot<-stack_assoc_plot(l, scores, cor,traits=paste0("Colocalisation between ", traits,". \n" , "cis-EQTL effect: ",finaldf_row$gene),type="prob",build=38,labels=top3,x.min=xmin,x.max=xmax)
    return(plot)
    }
  }
