


filter_genotype_matrix=function(genotypemat){

  tmp<-genotypemat

  message(paste0("Filtering genotype matrix. Starting with ",nrow(tmp)," snps."))

  tmp$counter_0<-rowSums(tmp==0)
  tmp$counter_1<-rowSums(tmp==1)
  tmp$counter_2<-rowSums(tmp==2)


  option1=tmp[tmp$counter_0>=2 & tmp$counter_1>=2,]
  option2=tmp[tmp$counter_0>=2 & tmp$counter_2>=2,]
  option3=tmp[tmp$counter_1>=2 & tmp$counter_2>=2,]

  intersected_rows<-unique(c(rownames(option1),rownames(option2),rownames(option3)))
  tmp<-tmp[intersected_rows,]

  message(paste0("Filtering complete. ",nrow(tmp), " snps have at least 2 individuals in at 2 out of the 3 dosage categories."))

  return(tmp)
}
