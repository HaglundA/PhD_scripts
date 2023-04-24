


library(ieugwasr)
library(parallel)

# make sure pvalue column is named "pval" and variant name is "rsid"



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
