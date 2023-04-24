


# Via seurat, mean counts

Idents(HIP_input_seuratObject) <- "Sample_ID"
HIP_count_ave <- AverageExpression(HIP_input_seuratObject, slot = "counts")
HIP_count_mean <- HIP_count_ave$RNA

# sum counts + cpm

### pseudobulk then CPM



# Input is list of filtered sc MX per cell type and metadata df
# Output is list of pseudobulk MX



psuedoBulk=function(scmatrix, metadata){

  CellTypes=names(scmatrix)
  pseudobulk=setNames(vector("list", length=length(CellTypes)),CellTypes)


  for (ct in CellTypes){
    pseudobulk[[ct]] <- Matrix.utils::aggregate.Matrix(t(scmatrix[[ct]]), groupings = metadata[colnames(scmatrix[[ct]]),"orig.ident"], fun = "sum")
    pseudobulk[[ct]]=log2(cpm(t(pseudobulk[[ct]]))+1)
  }

  saveRDS(pseudobulk, "pseudobulkSum_Log2CPM_filt1countPerSample.rds")

  return(pseudobulk)
}



agg_cpm=function(seuratlist){
  agg_count_list<-lapply(seuratlist,function(x){
    DefaultAssay(x)<-"RNA"
    metadata<-x[[]]

    counts<-as.data.frame(GetAssayData(x))
    agg<-Matrix.utils::aggregate.Matrix(t(counts), groupings = metadata[colnames(counts),"Sample_ID"], fun = "sum")
    return(t(agg))
  })

  # now normalise
  agg_count_list<-lapply(agg_count_list,function(x){
    norm_x<-log2(cpm(x)+1)
    return(norm_x)
  })

  names(agg_count_list)<-names(seuratlist)

  #now save matrices
  for(i in 1:length(agg_count_list)){
    write.table(agg_count_list[[i]],paste0(names(agg_count_list[i]),"_agg_cpm_pseudobulk.csv"))
  }

}
