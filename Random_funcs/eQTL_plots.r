


exp_mat<-read.table(exp_mat)
geno_mat<-read.table(geno_mat)





eQTL_plot_PD_covmat=function(exp_mat,geno_mat,snp,gene){


  covariates<-read.csv("/rds/general/user/ah3918/home/JAN21_WORK_DIR/Data/PD_Data/RAW.NONNORMALISED.COUNTS/MetaData.All_Samples1.csv",sep=",",row.names=1)

  tmp<-covariates[!duplicated(covariates$IndID),]
  tmp<-tmp[which(as.character(tmp$IndID) %in% colnames(exp_mat)),]
   tmp2<-data.frame(Gender=tmp$Gender,AgeAtDeath=tmp$AgeAtDeath,PMI=tmp$PMI)
   tmp2<-t(tmp2)
   colnames(tmp2)<-as.character(tmp$IndID)
   tmp2<-as.matrix(tmp2)
   covmat<-tmp2


   numericgeno<-as.matrix(geno_mat)
   numericexp<-as.matrix(exp_mat)

   lm_mat<-cbind(numericgeno[snp,],numericexp[gene,],covmat)
   colnames(lm_mat)<-c("geno","exp",colnames(covmat))
   lm_mat<-apply(lm_mat,2,as.numeric)


 res<-coef(summary(lm(geno~exp+Gender+PMI+AgeAtDeath, data=as.data.frame(lm_mat))))
 pval<-res["exp",4]
 pval<-round(pval,digits=6)
 return(pval)
}




get_gov_pval=function(exp_mat,geno_mat,snp,gene,covmat){

    common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
    geno_mat<-geno_mat[,common_names]
    exp_mat<-exp_mat[,common_names]
  covmat<-covmat[,common_names]
  covmat<-covmat[c("Age","PMI","Sample_Source","Sex"),]
  covmat<-t(covmat)
  #get LM p-value
  numericgeno<-as.matrix(geno_mat)
  numericexp<-as.matrix(exp_mat)

    lm_mat<-cbind(numericgeno[snp,],numericexp[gene,],covmat)
    colnames(lm_mat)<-c("geno","exp",colnames(covmat))
    lm_mat<-apply(lm_mat,2,as.numeric)


  res<-coef(summary(lm(geno~exp+Sex+PMI+Age+Sample_Source, data=as.data.frame(lm_mat))))
  pval<-res["exp",4]
  signif(pval, digits=3)

}

create_geno_plot=function(exp_mat,geno_mat,snp,gene,name){
  geno_mat=geno_mat[colnames(exp_mat)]


  #get LM p-value
  numericgeno<-as.matrix(geno_mat)
  numericexp<-as.matrix(exp_mat)
  sum<-summary(lm(numericgeno[snp,]~numericexp[gene,]))
  pval<-pf(sum$fstatistic[1],sum$fstatistic[2],sum$fstatistic[3],lower.tail=FALSE)
  attributes(pval)<-NULL
  pval<-round(pval,digits=4)

  df<-data.frame(snp=as.numeric(geno_mat[snp,]))
  df<-cbind(df,as.numeric(exp_mat[gene,]))
  colnames(df)<-c("geno","exp")

  test<-df$geno
  for(i in 1:length(test)){
      if(test[i]==0){
          test[i]<-"Homozygous_Alt"
      } else if(test[i]==1){
          test[i]<-"Heterozygous"
      }else{
          test[i]<-"Homozygous_Ref"
      }
  }
  df$geno<-test

  colnames(df)<-c("geno","exp")


  g<-ggplot(df,aes(y=exp,x=geno,fill=geno,group=geno)) + geom_boxplot(outlier.size=0, alpha=0.4, fill="grey")
  g<-g+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1))+annotate("text", x=3,size=5, y=0.01, label=paste0("p-val: ",pval))
  g<-g+ggtitle(paste0("eQTL effect of ",snp," on \n ",gene," in ",name))
  g<-g+theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x=element_text(size=12,face="bold"),axis.text.y=element_text(size=12,face="bold"))
  g+scale_x_discrete(limits=c("Homozygous_Ref","Heterozygous","Homozygous_Alt"))
}

add_pcs=
create_geno_cov_plot=function(exp_mat,
                              geno_mat,
                              snp,gene,
                              name,
                              covmat,
                            covs){

    common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
    geno_mat<-geno_mat[,common_names]
    exp_mat<-exp_mat[,common_names]

  # geno_mat=geno_mat[colnames(exp_mat)]
  covmat<-covmat[,colnames(exp_mat)]
  covmat<-covmat[covs,]
  covmat<-t(covmat)
  #get LM p-value
  numericgeno<-as.matrix(geno_mat)
  numericexp<-as.matrix(exp_mat)

  lm_mat<-cbind(numericgeno[snp,],numericexp[gene,],covmat)
  colnames(lm_mat)<-c("geno","exp",colnames(covmat))
  lm_mat<-apply(lm_mat,2,as.numeric)
  lm_model<-lm(paste("geno ~ exp+",paste(covs, collapse="+"),sep = ""), data=as.data.frame(lm_mat))

  res<-coef(summary(lm_model))
  pval<-res["exp",4]
    pval<-signif(pval, digits=3)
  resids<-resid(lm_model)

  df<-data.frame(snp=as.numeric(geno_mat[snp,]))
  df<-cbind(df,as.numeric(resids))
  colnames(df)<-c("geno","exp")

  test<-df$geno
  for(i in 1:length(test)){
      if(test[i]==0){
          test[i]<-"Homozygous_Alt"
      } else if(test[i]==1){
          test[i]<-"Heterozygous"
      }else{
          test[i]<-"Homozygous_Ref"
      }
  }
  df$geno<-test

  colnames(df)<-c("geno","exp")


  g<-ggplot(df,aes(y=exp,x=geno,fill=geno,group=geno)) + geom_boxplot(outlier.size=0, alpha=0.1, fill="grey")
  g<-g+geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5)+annotate("text", x=3,size=5, y=0.01, label=paste0("p-val: ",pval))
  g<-g+ggtitle(paste0("eQTL effect of ",snp," on \n ",gene," in ",name))
  g<-g+theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x=element_text(size=12,face="bold"),axis.text.y=element_text(size=12,face="bold"))
  g+scale_x_discrete(limits=c("Homozygous_Ref","Heterozygous","Homozygous_Alt"))
}



create_geno_cov_plot=function(exp_mat,
                              geno_mat,
                              snp,gene,
                              name,
                              covmat,
                            covs){

  #     common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
  #     geno_mat<-geno_mat[,common_names]
  #     exp_mat<-exp_mat[,common_names]

    # geno_mat=geno_mat[colnames(exp_mat)]

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
  # exp_model<-lm(paste("exp~",paste(covs, collapse="+"),sep = ""), data=as.data.frame(lm_mat))
   res<-coef(summary(lm_model))
    pval<-res["exp",4]
      pval<-signif(pval, digits=4)

    # adjusted_BF <- (coef(exp_model)[1] + resid(exp_model))

  #     test_mat<-as.data.frame(lm_mat)
  #     test_mat$adj_exp<-adjusted_BF

  #     lm_model<-lm(geno~adj_exp, data=test_mat)

  #   pval<-res["adj_exp",4]

    print(pval)
  resids<-resid(lm_model)


  df<-data.frame(snp=as.numeric(geno_mat[snp,]))
    df<-cbind(df,as.numeric(resids))
    colnames(df)<-c("geno","exp")

    test<-df$geno
    for(i in 1:length(test)){
        if(test[i]==0){
            test[i]<-"Alt"
        } else if(test[i]==1){
            test[i]<-"Het"
        }else{
            test[i]<-"Ref"
        }
    }
    df$geno<-test

    colnames(df)<-c("geno","exp")

    print(pval)
  g<-ggplot(df,aes(y=exp,x=geno,fill=geno,group=geno)) + geom_boxplot(outlier.size=0, alpha=0.1, fill="grey")
  g<-g+geom_dotplot(binaxis='y', stackdir='center',dotsize=0.8)+annotate("text", x=-Inf,size=5, y=Inf, label=paste0("p-val: ",pval))
  g<-g+ggtitle(paste0("eQTL effect of ",snp," on \n ",gene," in ",name))
  g<-g+theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x=element_text(size=12,face="bold"),axis.text.y=element_text(size=12,face="bold"))
  g+scale_x_discrete(limits=c("Ref","Het","Alt"))+scale_y_continuous(limits=c(-0.5,3))
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
  # exp_model<-lm(paste("exp~",paste(covs, collapse="+"),sep = ""), data=as.data.frame(lm_mat))


create_plots=function(lm_models)


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

        labs(subtitle = paste(" P =",signif(summary(fit)$coef[2,4], 3)))


    }

    plot_list<-list()
    plot_list[[1]]<-ggplotRegression(lm_model_list[[1]])+
    scale_fill_manual(values=cell_color_vectors[[1]])+
    theme_classic()+theme(legend.position="none")+
    ggtitle(cellnames[[1]])+
    scale_y_continuous(limits=c(lowerlimit,toplimit))+
    xlab("Genotype")+ylab("Expression")

    for(i in 2:length(cellnames)){
        plot_list[[i]]<-ggplotRegression(lm_model_list[[i]])+
    scale_fill_manual(values=cell_color_vectors[[i]])+
    theme_classic()+theme(legend.position="none",
                          axis.line.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          axis.text.y=element_blank(),
            axis.title.y=element_blank())+
        ggtitle(cellnames[[i]])+
        scale_y_continuous(limits=c(lowerlimit,toplimit))+
        xlab("Genotype")
        }
    plot_list[[5]]<-ggplotRegression(lm_model_list[[5]])+
    scale_fill_manual(values=cell_color_vectors[[5]])+
    theme_classic()+theme(legend.position="none")+
    ggtitle(cellnames[[5]])+
    scale_y_continuous(limits=c(lowerlimit,toplimit))+
    xlab("Genotype")+ylab("Expression")


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
  # exp_model<-lm(paste("exp~",paste(covs, collapse="+"),sep = ""), data=as.data.frame(lm_mat))
