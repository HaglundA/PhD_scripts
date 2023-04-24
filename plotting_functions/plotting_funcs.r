all_plots_violin=function(exp_mat_list,exp_loc,geno_mat,covmat,gene,snp,covs,allele_mat){



    lm_model_list<-list()

    #limits for plots
    toplimit<-vector()
    lowerlimit<-vector()
    cellnames<-names(exp_mat_list)
    for(i in 1:length(exp_mat_list)){

    # cellnames<-names(exp_mat_list[i])

    exp_mat<-exp_mat_list[[i]]
    exp_mat<-filter_pseudobulk(exp_mat,minimum_indivs = 3)
    exp_mat<-exp_mat[rownames(exp_mat) %in% exp_loc$geneid,]
    common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
    geno_mat_tmp<-geno_mat[common_names]
    exp_mat<-exp_mat[common_names]

    covmat_tmp<-covmat[,colnames(exp_mat)]
    covmat_tmp<-add_pcs(covmat = covmat_tmp,exp_mat=exp_mat,npcs = 51)

    snp=snp
    gene=gene
    covs=covs
    lm_model_list[[i]]<-return_lm_model_covs(exp_mat,geno_mat_tmp,snp,gene,covmat=covmat_tmp,covs=covs)

    toplimit<-c(toplimit,max(lm_model_list[[i]]$model$exp))
       lowerlimit<-c(lowerlimit,min(lm_model_list[[i]]$model$exp))
    }

    #plot prep work
    #set limits
    toplimit<-max(toplimit)
    toplimit<-toplimit+toplimit*0.05
    lowerlimit<-min(lowerlimit)
    lowerlimit<-lowerlimit-lowerlimit*0.05

    #get alleles
    allele_mat<-allele_mat[which(allele_mat$snp==snp),]
    ref<-paste0(allele_mat$ref[1],"/",allele_mat$ref[1])
    alt<-paste0(allele_mat$alt[1],"/",allele_mat$alt[1])
    het<-paste0(allele_mat$alt[1],"/",allele_mat$ref[1])


    ### plotting

    ggplotRegression <- function (fit) {



    require(ggplot2)
    test<-fit$model$geno
          for(i in 1:length(test)){
          if(test[i]==0){
              test[i]<-alt
          } else if(test[i]==1){
              test[i]<-het
          }else{
              test[i]<-ref
          }
      }
        fit$model$geno_char<-test
        df<-fit$model

        smooth<-predict(fit)

    ggplot() +
        geom_violin(df, mapping=aes(x = geno_char, y = exp,fill=geno_char),alpha=0.9)+
    geom_point(df, mapping=aes(x = geno_char, y = exp,fill=geno_char))+
        stat_smooth(method="lm",col="red",aes(x = as.numeric(as.factor(geno_char)), y = exp),data=df,se=F)+

        labs(subtitle = paste(" P =",signif(summary(fit)$coef[2,4], 5)))


    }

    plot_list<-list()
    plot_list[[1]]<-ggplotRegression(lm_model_list[[1]])+
    scale_fill_manual(values=cell_color_vectors[[1]])+
    theme_classic()+theme(legend.position="none")+
    ggtitle(cellnames[[1]])+scale_y_continuous(limits=c(lowerlimit,toplimit))+xlab("Genotype")+ylab("Expression")

    for(i in 2:length(cellnames)){
        plot_list[[i]]<-ggplotRegression(lm_model_list[[i]])+
    scale_fill_manual(values=cell_color_vectors[[i]])+
    theme_classic()+theme(legend.position="none",axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.line.y=element_blank())+
        ggtitle(cellnames[[i]])+scale_y_continuous(limits=c(lowerlimit,toplimit))+xlab("Genotype")
        }
    plot_list[[5]]<-ggplotRegression(lm_model_list[[5]])+
    scale_fill_manual(values=cell_color_vectors[[5]])+
    theme_classic()+theme(legend.position="none")+
    ggtitle(cellnames[[5]])+scale_y_continuous(limits=c(lowerlimit,toplimit))+xlab("Genotype")+ylab("Expression")


    ##return final plot
  title <- ggdraw() + draw_label(paste0("Effect of ",snp," on ",gene), fontface='bold',size=17)
    grid1<-suppressMessages(plot_grid(plot_list[[1]],
                              plot_list[[2]],
                              plot_list[[3]],
                              plot_list[[4]],
                              plot_list[[5]],
                              plot_list[[6]],plot_list[[7]],plot_list[[8]],nrow=2,align="h"))

    plot_grid(title,grid1,align="h",ncol=1,rel_heights =c(0.2,2) )


    }
all_plots_boxplot=function(exp_mat_list,exp_loc,geno_mat,covmat,gene,snp,covs,allele_mat){



    lm_model_list<-list()

    #limits for plots
    toplimit<-vector()
    lowerlimit<-vector()
    cellnames<-names(exp_mat_list)

    if(any(colnames(geno_mat)=="S03.009.B")==TRUE){
      names(geno_mat)[names(geno_mat)=="S03.009.B"]<-"S03.009"
    }
    for(i in 1:length(exp_mat_list)){

      celltype<-names(exp_mat_list[i])

        exp_mat<-exp_mat_list[[i]]
        exp_mat<-filter_pseudobulk(exp_mat,minimum_indivs = 3)
        exp_mat<-exp_mat[rownames(exp_mat) %in% exp_loc$geneid,]
        common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
        geno_mat_tmp<-geno_mat[common_names]
        exp_mat<-exp_mat[common_names]

        covmat_tmp<-covmat[,colnames(exp_mat)]
        covmat_tmp<-add_pcs(covmat = covmat_tmp,exp_mat=exp_mat,npcs = 51)
    snp=snp
    gene=gene
    covs=covs
    lm_model_list[[i]]<-return_lm_model_covs(exp_mat,geno_mat_tmp,snp,gene,covmat=covmat_tmp,covs=covs)

    toplimit<-c(toplimit,max(lm_model_list[[i]]$model$exp))
       lowerlimit<-c(lowerlimit,min(lm_model_list[[i]]$model$exp))
    }

    #plot prep work
    #set limits
    toplimit<-max(toplimit)
    toplimit<-toplimit+toplimit*0.05
    lowerlimit<-min(lowerlimit)
    lowerlimit<-lowerlimit-lowerlimit*0.05

    #get alleles
    allele_mat<-allele_mat[which(allele_mat$snp==snp),]
    ref<-paste0(allele_mat$ref[1],"/",allele_mat$ref[1])
    alt<-paste0(allele_mat$alt[1],"/",allele_mat$alt[1])
    het<-paste0(allele_mat$alt[1],"/",allele_mat$ref[1])



    ##colors

    darker.col = function(color, how.much = 30){
      colorRampPalette(c(color, "black"))(100)[how.much]
    }
    lighter.col = function(color, how.much = 30){
      colorRampPalette(c(color, "white"))(100)[how.much]
    }
    cell_color_vectors<-list()

    for(i in 1:length(cellnames)){
        color=ggsci::pal_nejm("default")(8)[i]
        cell_color_vectors[[i]]<-rev(c(color,lighter.col(color,20),lighter.col(color,50)))
        }


    ### plotting

    ggplotRegression <- function (fit) {



    require(ggplot2)
    test<-fit$model$geno
          for(i in 1:length(test)){
          if(test[i]==0){
              test[i]<-alt
          } else if(test[i]==1){
              test[i]<-het
          }else{
              test[i]<-ref
          }
      }
        fit$model$geno_char<-test
        df<-fit$model

        smooth<-predict(fit)

    ggplot() +
        geom_boxplot(df, mapping=aes(x = geno_char, y = exp,fill=geno_char),alpha=0.9)+
    geom_point(df, mapping=aes(x = geno_char, y = exp,fill=geno_char))+
        stat_smooth(method="lm",col="red",aes(x = as.numeric(as.factor(geno_char)), y = exp),data=df,se=F)+

        labs(subtitle = paste(" P =",signif(summary(fit)$coef[2,4], 5)))


    }

    plot_list<-list()
    plot_list[[1]]<-ggplotRegression(lm_model_list[[1]])+
    scale_fill_manual(values=cell_color_vectors[[1]])+
    theme_classic()+theme(legend.position="none")+
    ggtitle(cellnames[[1]])+scale_y_continuous(limits=c(lowerlimit,toplimit))+xlab("Genotype")+ylab("Expression")

    for(i in 2:length(cellnames)){
        plot_list[[i]]<-ggplotRegression(lm_model_list[[i]])+
    scale_fill_manual(values=cell_color_vectors[[i]])+
    theme_classic()+theme(legend.position="none",axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.line.y=element_blank())+
        ggtitle(cellnames[[i]])+scale_y_continuous(limits=c(lowerlimit,toplimit))+xlab("Genotype")
        }
    plot_list[[5]]<-ggplotRegression(lm_model_list[[5]])+
    scale_fill_manual(values=cell_color_vectors[[5]])+
    theme_classic()+theme(legend.position="none")+
    ggtitle(cellnames[[5]])+scale_y_continuous(limits=c(lowerlimit,toplimit))+xlab("Genotype")+ylab("Expression")


    ##return final plot
  title <- cowplot::ggdraw() + draw_label(paste0("Effect of ",snp," on ",gene), fontface='bold',size=17)
    grid1<-suppressMessages(plot_grid(plot_list[[1]],
                              plot_list[[2]],
                              plot_list[[3]],
                              plot_list[[4]],
                              plot_list[[5]],
                              plot_list[[6]],plot_list[[7]],plot_list[[8]],nrow=2,align="h"))

    plot_grid(title,grid1,align="h",ncol=1,rel_heights =c(0.2,2) )

}


filter_pseudobulk=function(exp_mat,minimum_indivs){


  exp_mat<-exp_mat[rowSums(exp_mat>0)>=minimum_indivs,]
  return(exp_mat)



}
add_pcs=function(covmat,exp_mat,npcs){
    covs<-covmat    
    #     covs<-covs[,colnames(exp_mat)]
      
      pcs<-prcomp(exp_mat,scale=T,center=T)
      pcs<-pcs$rotation
      pcs<-pcs[,1:npcs]
      pcs<-t(pcs)
      pcs<-pcs[,colnames(covs)]
      rownames(pcs)<-paste0(rep("PC."),1:npcs)
      covs<-as.matrix(covs)
      covs<-rbind(covs,pcs)
    return(covs)


}
return_lm_model_covs=function(exp_mat,





                              geno_mat,
                              snp,gene,
                              name,
                              covmat,
                            covs){

    covmat<-covmat[covs,]
  #   covmat<-covmat[,colnames(exp_mat)]
    covmat<-t(covmat)
    #get LM p-value
    numericgeno<-as.matrix(geno_mat)
    numericexp<-as.matrix(exp_mat)

    lm_mat<-cbind(numericgeno[snp,],numericexp[gene,],covmat)
    colnames(lm_mat)<-c("geno","exp",colnames(covmat))
    lm_mat<-apply(lm_mat,2,as.numeric)

    lm_model<-lm(paste("geno~exp +",paste(covs, collapse="+"),sep = ""), data=as.data.frame(lm_mat))
    return(lm_model)
}


 