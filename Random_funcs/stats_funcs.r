




eQTL_lm_pval=function(exp_mat,geno_mat,snp,gene,cov_mat){

    cov_mat<-cov_mat[colnames(exp_mat)]
  cov_mat<-t(cov_mat)
  geno_mat<-geno_mat[colnames(exp_mat)]


  numericgeno<-as.matrix(geno_mat)
  numericexp<-as.matrix(exp_mat)
  lm_mat<-cbind(numericgeno[snp,],numericexp[gene,],cov_mat)
  colnames(lm_mat)<-c("geno","exp",colnames(cov_mat))

  res<-coef(summary(lm(geno~exp+AgeAtDeath+Gender+PMI, data=as.data.frame(lm_mat))))
  pval<-res["exp",4]
  pval<-round(pval,digits=6)
  return(pval)


}


eQTL_rlm_pval=function(exp_mat,geno_mat,snp,gene,cov_mat){
   cov_mat<-cov_mat[colnames(exp_mat)]
  cov_mat<-t(cov_mat)
  geno_mat<-geno_mat[colnames(exp_mat)]


  numericgeno<-as.matrix(geno_mat)
  numericexp<-as.matrix(exp_mat)

  lm_mat<-cbind(numericgeno[snp,],numericexp[gene,],cov_mat)
  colnames(lm_mat)<-c("geno","exp",colnames(cov_mat))

  res<-rlm(geno~exp+AgeAtDeath+Gender+PMI, data=as.data.frame(lm_mat))
  pval<-f.robftest(res,var="exp")$p.value
  return(pval)
}
