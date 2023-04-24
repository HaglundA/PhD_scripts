


filter_pseudobulk=function(exp_mat,minimum_indivs){


  exp_mat<-exp_mat[rowSums(exp_mat>0)>=minimum_indivs,]
  return(exp_mat)



}
