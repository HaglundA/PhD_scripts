




MatrixEQTL_out_list<-list.files("../",pattern="MatrixEQTL",full.names=T)
mateqtlouts<-lapply(MatrixEQTL_out_list,read.table,head=T)
cellnames<-gsub("_TreeQTL_eAssoc_cis_by_snp.txt","",list.files("../",pattern="TreeQTL"))





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

    fdr_0_0_5<-x[x$FDR<0.05,]

    # n_snp_genepairs_fdr_0_0_5<-nrow(fdr_0_0_5)
    n_eGenes_fdr_0_0_5<-length(unique(fdr_0_0_5$gene))
    n_eSNPS_fdr_0_0_5<-length(unique(fdr_0_0_5$SNP))

    tmp<-c(n_eSNPs,
      n_eGenes,
      n_eSNPS_fdr_0_2,
      n_eGenes_fdr_0_2,
      n_eSNPS_fdr_0_0_5,
        n_eGenes_fdr_0_0_5)

        names(tmp)<-c("n_eSNPs",
        "n_eGenes",
        "n_eSNPs_fdr_0_2","n_eGenes_fdr_0_2","n_eSNPs_fdr_0_0_5","n_eGenes_fdr_0_0_5")

        return(tmp)


}

mateqtls_filt<-lapply(mateqtlouts,function(x){
    x<-x[x$p.value<0.05,]
    return(x)
})

allres<-lapply(mateqtls_filt,get_info)

test<-do.call(cbind,allres)
options(repr.plot.width=NULL,repr.plot.height=NULL)
my_bar<-barplot(test[4:6,],col=c("grey24","grey55","grey75"),main="All cis-eQTLs, FDR < 20%", las=2)
legend("topleft",c("eGenes","eSNPs","SNP/Gene pairs"),col=c("grey75","grey55","grey24"),lwd=10)
