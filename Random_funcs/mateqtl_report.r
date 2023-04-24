


summarize_outputs=function(outdir){

  setwd(outdir)
  cellnames<-gsub("_cis_MatrixEQTLout.txt","",list.files())
  message("Reading in MatrixEQTL outputs..")
  mateqtlouts<-lapply(list.files(pattern="_cis_MatrixEQTLout"),read.table)
  names(mateqtlouts)<-cellnames

  message("Saving MatrixEQTL outputs as list object..")
  saveRDS(mateqtlouts,"mateqtlouts.rds")

  res<-lapply(mateqtlouts,get_info)
  res<-as.data.frame(do.call(cbind,res))
  colnames(res)<-cellnames

  unique<-get_unique(mateqtlouts)
  res$total<-rowSums(res)
  res$unique<-unique

  write.table(res,"summary_eQTL_numbers.txt"))


}
get_info=function(x){
  tmp<-vector()


    n_snp_genepairs<-nrow(x)
    n_eGenes<-length(unique(x$gene))
    n_eSNPs<-length(unique(x$SNP))

    # same again for FDR
    fdr_0_2<-x[x$FDR<0.2,]

    # n_snp_genepairs_fdr_0_2<-nrow(fdr_0_2)
    n_eGenes_fdr_0_2<-length(unique(fdr_0_2$gene))
    n_eSNPS_fdr_0_2<-length(unique(fdr_0_2$SNP))

    fdr_0_1<-x[x$FDR<0.1,]

    # n_snp_genepairs_fdr_0_0_5<-nrow(fdr_0_0_5)
    n_eGenes_fdr_0_1<-length(unique(fdr_0_1$gene))
    n_eSNPS_fdr_0_1<-length(unique(fdr_0_1$SNP))

    fdr_0_0_5<-x[x$FDR<0.05,]

    # n_snp_genepairs_fdr_0_0_5<-nrow(fdr_0_0_5)
    n_eGenes_fdr_0_0_5<-length(unique(fdr_0_0_5$gene))
    n_eSNPS_fdr_0_0_5<-length(unique(fdr_0_0_5$SNP))

    tmp<-c(n_eSNPs,
      n_eGenes,
      n_eSNPS_fdr_0_2,
      n_eGenes_fdr_0_2,
      n_eSNPS_fdr_0_1,
      n_eGenes_fdr_0_1,
      n_eSNPS_fdr_0_0_5,
        n_eGenes_fdr_0_0_5)

        names(tmp)<-c("n_eSNPs",
        "n_eGenes",
        "n_eSNPs_fdr_0_2","n_eGenes_fdr_0_2",
        "n_eSNPS_fdr_0_1","n_eSNPS_fdr_0_1",
        "n_eSNPs_fdr_0_0_5","n_eGenes_fdr_0_0_5")

        return(tmp)


}

get_unique=function(x){

    #x is a mateqtlouts list

    genes<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes<-length(unique(unlist(genes)))
    snps<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps<-length(unique(unlist(snps)))

    #20%fdr------

    x<-lapply(x,function(a){
        a<-a[a$FDR<0.2,]
    })
    genes_0_2_FDR<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes_0_2_FDR<-length(unique(unlist(genes_0_2_FDR)))
    snps_0_2_FDR<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps_0_2_FDR<-length(unique(unlist(snps_0_2_FDR)))

    #10%fdr------


    x<-lapply(x,function(a){
        a<-a[a$FDR<0.1,]
    })
    genes_0_1_FDR<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes_0_1_FDR<-length(unique(unlist(genes_0_1_FDR)))
    snps_0_1_FDR<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps_0_1_FDR<-length(unique(unlist(snps_0_1_FDR)))




    #5%fdr------

    x<-lapply(x,function(a){
        a<-a[a$FDR<0.05,]
    })
    genes_0_0_5_FDR<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes_0_0_5_FDR<-length(unique(unlist(genes_0_0_5_FDR)))
     snps_0_0_5_FDR<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps_0_0_5_FDR<-length(unique(unlist(snps_0_0_5_FDR)))

    finalvec<-c(snps,genes,snps_0_2_FDR,genes_0_2_FDR,snps_0_1_FDR,genes_0_1_FDR,
      snps_0_0_5_FDR,genes_0_0_5_FDR)
    return(finalvec)


}
