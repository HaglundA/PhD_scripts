





compare_rank=function(bryois_all,mateqtlouts,celltype){

  Oligo<-bryois_all
  Oligo<-Oligo[Oligo$pvalue<0.05,]
  Oligo<-Oligo[order(Oligo$p.value),]
  colnames(Oligo)<-c("gene_ensg","snp","direction","pvalue","beta")
  Oligo$gene<-do.call(rbind,strsplit(Oligo$gene_ensg,"_"))[,1]

  tmp<-mateqtlouts[[celltype]]
  tmp<-tmp[tmp$p.value<0.05,]
  tmp<-tmp[order(tmp$p.value),]

  tmp<-tmp[match(Oligo$gene,tmp$gene),]
  tmp<-tmp[match(Oligo$snp)]



}




compare_all_snpgenes=function(biglist,
  celltype,
  # bryois_df,
  psychencode,
  gtex_pfc,
  gtex_hip,
  metabrain){

    #first, Bryois
    # astro<-bryois_df
    # astro$V1<-as.character(astro$V1)
    # astro$V2<-as.character(astro$V2)
    # bryois_astro<-astro
    # colnames(bryois_astro)<-c("gene","SNP","V3","p.value","beta")
    # genes<-bryois_astro$gene
    # genes<-strsplit(genes,"_")
    # genes<-as.data.frame(do.call(rbind,genes))
    # bryois_astro$gene<-genes$V1
    # bryois_astro$SNP_gene<-paste0(bryois_astro$gene,"_",bryois_astro$SNP)
    # bryois_astro<-bryois_astro[order(bryois_astro$p.value),]
    # tmp_bryois<-bryois_astro[bryois_astro$p.value<0.05,]
    #
    #
    intersected<-vector()
    total<-vector()
    # for(i in 1:length(biglist)){
    #     tmp<-biglist[[i]]
    #     tmp<-tmp[[celltype]]
    #     tmp<-tmp[tmp$FDR<0.05,]
    #     tmp$SNP_gene<-paste0(tmp$gene,"_",tmp$SNP)
    #
    #     total<-c(total,nrow(tmp))
    #
    #     intersected<-c(intersected,length(intersect(tmp$SNP_gene,tmp_bryois$SNP_gene)))
    #   }


    ## Psychencode
    test<-psychencode
    for(i in 1:length(biglist)){
        tmp<-biglist[[i]]
        tmp<-tmp[[celltype]]
        tmp<-tmp[tmp$FDR<0.05,]
        tmp$SNP_gene<-paste0(tmp$gene,"_",tmp$SNP)


        total<-c(total,nrow(tmp))
        intersected<-c(intersected,length(intersect(tmp$SNP_gene,test$SNP_gene)))


    }


    # gtex_pfc

    for(i in 1:length(biglist)){
        tmp<-biglist[[i]]
        tmp<-tmp[[celltype]]
        tmp<-tmp[tmp$FDR<0.05,]
        tmp$SNP_gene<-paste0(tmp$gene,"_",tmp$SNP)



        intersected<-c(intersected,length(intersect(tmp$SNP_gene,gtex_pfc$SNP_gene)))


    }

    # gtex_pfc

    for(i in 1:length(biglist)){
        tmp<-biglist[[i]]
        tmp<-tmp[[celltype]]
        tmp<-tmp[tmp$FDR<0.05,]
        tmp$SNP_gene<-paste0(tmp$gene,"_",tmp$SNP)



        intersected<-c(intersected,length(intersect(tmp$SNP_gene,gtex_hip$SNP_gene)))


    }

    #metabrain

    for(i in 1:length(biglist)){
        tmp<-biglist[[i]]
        tmp<-tmp[[celltype]]
        tmp<-tmp[tmp$FDR<0.05,]
        tmp$SNP_gene<-paste0(tmp$gene,"_",tmp$SNP)



        intersected<-c(intersected,length(intersect(tmp$SNP_gene,metabrain$SNP_gene)))


    }




    pcs<-c(0,10,20,30,40,50)
    percent<-round(intersected/total,digits = 3)
    percent<-c(total,percent)
    master_vec<-c(total,intersected)
    category<-c(rep("total",length(biglist)),
    rep("Psychencode",length(biglist)),rep("GTEx_cortex",length(biglist)),
  rep("GTEx_hip",length(biglist)),rep("Metabrain_cortex",length(biglist)))
    return(list(pcs,master_vec,category,percent))







}






compare_all_egenes=function(biglist,
    celltype,
    # bryois_df,
    psychencode,
    gtex_pfc,
    gtex_hip,
    metabrain){

      #first, Bryois
      # astro<-bryois_df
      # astro$V1<-as.character(astro$V1)
      # astro$V2<-as.character(astro$V2)
      # bryois_astro<-astro
      # colnames(bryois_astro)<-c("gene","SNP","V3","p.value","beta")
      # genes<-bryois_astro$gene
      # genes<-strsplit(genes,"_")
      # genes<-as.data.frame(do.call(rbind,genes))
      # bryois_astro$gene<-genes$V1
      # bryois_astro$gene<-paste0(bryois_astro$gene,"_",bryois_astro$SNP)
      # bryois_astro<-bryois_astro[order(bryois_astro$p.value),]
      # tmp_bryois<-bryois_astro[bryois_astro$p.value<0.05,]
      #
      #
      intersected<-vector()
      total<-vector()
      # for(i in 1:length(biglist)){
      #     tmp<-biglist[[i]]
      #     tmp<-tmp[[celltype]]
      #     tmp<-tmp[tmp$FDR<0.05,]
      #     tmp$gene<-paste0(tmp$gene,"_",tmp$SNP)
      #
      #     total<-c(total,nrow(tmp))
      #
      #     intersected<-c(intersected,length(intersect(tmp$gene,tmp_bryois$gene)))
      #   }


      ## Psychencode
      test<-psychencode
      test$gene<-as.character(test$gene)
      for(i in 1:length(biglist)){
          tmp<-biglist[[i]]
          tmp<-tmp[[celltype]]
          tmp<-tmp[tmp$FDR<0.05,]
          # tmp$gene<-paste0(tmp$gene,"_",tmp$SNP)


          total<-c(total,length(unique(tmp$gene)))
          intersected<-c(intersected,length(intersect(unique(tmp$gene),unique(test$gene))))


      }


      # gtex_pfc
      gtex_pfc$gene_name<-as.character(gtex_pfc$gene_name)
      for(i in 1:length(biglist)){
          tmp<-biglist[[i]]
          tmp<-tmp[[celltype]]
          tmp<-tmp[tmp$FDR<0.05,]
          # tmp$gene<-paste0(tmp$gene,"_",tmp$SNP)



          intersected<-c(intersected,length(intersect(unique(tmp$gene),unique(gtex_pfc$gene_name))))


      }

      # gtex_pfc
      gtex_hip$gene_name<-as.character(gtex_hip$gene_name)
      for(i in 1:length(biglist)){
          tmp<-biglist[[i]]
          tmp<-tmp[[celltype]]
          tmp<-tmp[tmp$FDR<0.05,]
          # tmp$gene<-paste0(tmp$gene,"_",tmp$SNP)



          intersected<-c(intersected,length(intersect(unique(tmp$gene),unique(gtex_hip$gene_name))))


      }

      #metabrain

      for(i in 1:length(biglist)){
          tmp<-biglist[[i]]
          tmp<-tmp[[celltype]]
          tmp<-tmp[tmp$FDR<0.05,]
          tmp$SNP_gene<-paste0(tmp$gene,"_",tmp$SNP)



          intersected<-c(intersected,length(intersect(unique(tmp$gene),unique(metabrain$Gene.Symbol))))


      }



      pcs<-c(0,10,20,30,40,50)
      percent<-round(intersected/total,digits = 3)
      percent<-c(total,percent)
      master_vec<-c(total,intersected)
      category<-c(rep("total",length(biglist)),
      rep("Psychencode",length(biglist)),rep("GTEx_cortex",length(biglist)),
    rep("GTEx_hip",length(biglist)),rep("Metabrain_cortex",length(biglist)))
      return(list(pcs,master_vec,category,percent))







}
