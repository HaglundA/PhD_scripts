


treeqtls<-lapply(list.files("../",pattern="TreeQTL",full.names=T),read.table,head=T)
cellnames<-gsub("_TreeQTL_eAssoc_cis_by_snp.txt","",list.files("../",pattern="TreeQTL"))
names(treeqtls)<-cellnames


create_upsetR=function(treeqtl_list){
  treeqtls<-lapply(treeqtls,function(x){
    x$SNP_genepair<-paste0(x$SNP,"_",x$gene)
    return(x)
  })

  SNP_genepairlist<-lapply(treeqtls,function(x){
    return(x$SNP_genepair)
  })
  genelist<-lapply(treeqtls,function(x){
    return(unique(x$gene))
  })

  names(SNP_genepairlist)<-names(treeqtls)
  names(genelist)<-names(treeqtls)
  upset(fromList(SNP_genepair,nsets=length(treeqtls)))
  upset(fromList(SNP_genepairlist),nsets=length(treeqtls)))
}



create_upsetR=function(mateqtlouts){

  treeqtls<-mateqtlouts

  SNP_genepairlist<-lapply(treeqtls,function(x){
    return(rownames(x))
  })
  genelist<-lapply(treeqtls,function(x){
    return(unique(x$gene))
  })

  names(SNP_genepairlist)<-names(treeqtls)
  names(genelist)<-names(treeqtls)
  upset(fromList(SNP_genepair,nsets=length(treeqtls)))
  upset(fromList(SNP_genepairlist),nsets=length(treeqtls)))
}
