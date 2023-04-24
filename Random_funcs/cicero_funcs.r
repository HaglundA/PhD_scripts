





# get peak from atac matrix given a snp or gene


get_peak=function(chr,pos,atac_mat){

  test<-rownames(atac_mat)
  test<-strsplit(test,"-")
  test<-do.call(rbind,test)
  colnames(test)<-c("chrom","start","end")
  test<-as.data.frame(test)

  test$start<-as.numeric(test$start)
  test$end<-as.numeric(test$end)
  test<-test[which(test$chr==chr),]

  test<-subset(test,start<pos & end>pos)

  return(paste0(test$chrom,"-",test$start,"-",test$end))

}






get_tss_pos_gene=function(gene_anno,gene){

  gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
  gene_anno$gene <- gene_anno$gene_id
  gene_anno$transcript <- gene_anno$transcript_id
  gene_anno$symbol <- gene_anno$gene_name

  gene_anno<-gene_anno[,c("chromosome","start","end","strand","type","gene","transcript","symbol")]
  pos <- subset(gene_anno, strand == "+")
  pos <- pos[order(pos$start),]
  pos <- pos[!duplicated(pos$transcript),] # remove all but the first exons per transcript
  pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS

  neg <- subset(gene_anno, strand == "-")
  neg <- neg[order(neg$start, decreasing = TRUE),]
  neg <- neg[!duplicated(neg$transcript),] # remove all but the first exons per transcript
  neg$start <- neg$end - 1

  gene_annotation_sub <- rbind(pos, neg)
  gene_annotation_sub <- gene_annotation_sub[,c(1:3, 8)]
  names(gene_annotation_sub)[4] <- "gene"
  colnames(gene_annotation_sub)<-c("chrom","position","end","gene")

  # overlap the gene
  gene_df<-gene_annotation_sub[which(gene_annotation_sub$gene==gene),]

 return(gene_df)

}



get_atac_correlation=function(snp,snplocs,gene,atac_mat,gene_anno){


  snplocs<-snplocs[which(snplocs$annot==snp),]

  snp_chr<-snplocs$chrom
  snp_pos<-snplocs$position

  snp_peak<-get_peak(snp_chr,snp_pos,atac_mat)

  gene_tss_pos<-get_tss_pos_gene(gene_anno,gene)

  gene_chr<-gene_tss_pos$chrom[1]
  gene_pos<-gene_tss_pos$position[1]

  gene_peak<-get_peak(gene_chr,gene_pos,atac_mat)

  tmp<-rbind(gene_peak,snp_peak)
  return(tmp)
}



create_atac_cds=function(indata){

  fixednames<-gsub("-","_",rownames(indata))
  rownames(indata)<-fixednames

  cellinfo<-as.data.frame(colnames(indata))
  names(cellinfo) <- "cells"
  row.names(cellinfo) <- cellinfo$cells


  peakinfo<-rownames(indata)
  df<-do.call(rbind,strsplit(peakinfo,"_"))
  df<-as.data.frame(df)
  rownames(df)<-peakinfo
  df$site_name<-peakinfo
  peakinfo<-df
  colnames(peakinfo)<-c("chr","bp1","bp2","site_name")

  fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
  pd <- methods::new("AnnotatedDataFrame", data = cellinfo)


  input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
  return(input_cds)

}
create_exp_cds=function(inputmat){

  exp_cells<-data.frame(cells=colnames(inputmat))
  rownames(exp_cells)<-colnames(inputmat)
  head(exp_cells)

  exp_info<-data.frame(gene=rownames(inputmat),gene_2=rownames(inputmat))
  rownames(exp_info)<-rownames(inputmat)
  head(exp_info)

  fd <- methods::new("AnnotatedDataFrame", data = exp_info)
  pd <- methods::new("AnnotatedDataFrame", data = exp_cells)


   exp_cds <-  suppressWarnings(newCellDataSet(inputmat,
                             phenoData = pd,
                             featureData = fd,
                             expressionFamily=VGAM::binomialff(),
                             lowerDetectionLimit=0))
   return(exp_cds)

}

get_atac_coords=function(atac_cds){
  input_cds<-atac_cds
  input_cds@expressionFamily@vfamily <- "binomialff"
  input_cds <- monocle::detectGenes(input_cds)

  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]


  set.seed(2017)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                      reduction_method = 'tSNE', norm_method = "none",check_duplicates=FALSE)
  tsne_coords <- t(reducedDimA(input_cds))
  row.names(tsne_coords) <- row.names(pData(input_cds))
  return(tsne_coords)
}

create_bagged_pseuduobulk=function(cds,
  tsne_dims,
  k_bags=50,
  summary_stats=NULL,
  size_factor_normalize=TRUE,
  seed=2016,
  silent=FALSE){

  set.seed(seed)
  reduced_coordinates<-tsne_dims
  k<-k_bags
  reduced_coordinates <- as.data.frame(reduced_coordinates) #dim 1
  reduced_coordinates <- reduced_coordinates[colnames(cds),]# Create a k-nearest neighbors map
   nn_map <- FNN::knn.index(reduced_coordinates, k=(k-1)) # no data.frame wrapper

   # nn_map dims<-1014 - 49

   # bind cellnames and make data frame
   row.names(nn_map) <- row.names(reduced_coordinates)
   nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))

   # define number of choices
   good_choices <- seq_len(nrow(nn_map))

   # sample first
   choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
   chosen <- good_choices[choice]
   good_choices <- good_choices[good_choices != good_choices[choice]]
   it <- 0
   k2 <- k * 2 # Compute once

   # function for sapply
   get_shared <- function(other, this_choice) {
     k2 - length(union(cell_sample[other,], this_choice))
   }

   while (length(good_choices) > 0 & it < 5000) { # slow
     it <- it + 1
     choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
     new_chosen <- c(chosen, good_choices[choice])
     good_choices <- good_choices[good_choices != good_choices[choice]]
     cell_sample <- nn_map[new_chosen,]

     others <- seq_len(nrow(cell_sample) - 1)
     this_choice <- cell_sample[nrow(cell_sample),]
     shared <- sapply(others, get_shared, this_choice = this_choice)

     if (max(shared) < .9 * k) {
       chosen <- new_chosen
     }
   }

   cell_sample <- nn_map[chosen,]

   if(!silent) {
     # Only need this slow step if !silent
     combs <- combn(nrow(cell_sample), 2)

     shared <- apply(combs, 2, function(x) {  #slow
       k2 - length(unique(as.vector(cell_sample[x,])))
     })

     # message(paste0("Overlap QC metrics:\nCells per bin: ", k,
     #                "\nMaximum shared cells bin-bin: ", max(shared),
     #                "\nMean shared cells bin-bin: ", mean(shared),
     #                "\nMedian shared cells bin-bin: ", median(shared)))

     if (mean(shared)/k > .1) warning("On average, more than 10% of cells are shared between paired bins.")
   }

   exprs_old <- exprs(cds)

   mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in% cell_sample[x,,drop=FALSE])
   mask <- Matrix::Matrix(mask)
   new_exprs <- exprs_old %*% mask

   new_exprs <- Matrix::t(new_exprs)
   new_exprs <- as.matrix(new_exprs)
   pdata <- pData(cds)
   new_pcols <- "agg_cell"
   if(!is.null(summary_stats)) {
     new_pcols <- c(new_pcols, paste0("mean_",summary_stats))
   }

   new_pdata <- plyr::adply(cell_sample,1, function(x) {
     sub <- pdata[x,]
     df_l <- list()
     df_l["temp"] <- 1
     for (att in summary_stats) {
       df_l[paste0("mean_", att)] <- mean(sub[,att])
     }
     data.frame(df_l)
   })

   new_pdata$agg_cell <- paste("agg", chosen, sep="")
   new_pdata <- new_pdata[,new_pcols, drop = FALSE] # fixes order, drops X1 and temp

   row.names(new_pdata) <- new_pdata$agg_cell
   row.names(new_exprs) <- new_pdata$agg_cell
   new_exprs <- as.matrix(t(new_exprs))



  return(new_exprs)
}


#cicero source









generate_windows <- function(window, genomic_coords) {
  if(!is(genomic_coords, "data.frame")) {
    chr_maxes <- read.table(genomic_coords)
  } else {
    chr_maxes <- genomic_coords
  }
  names(chr_maxes) <- c("V1", "V2")
  win_ranges <- plyr::ddply(chr_maxes, plyr::.(V1), function(x) {
    r <- seq(from = 1, to = x$V2[1], by = window/2)
    l <- r + window - 1
    data.frame(start = r, end = l)
  })
  gr <- GenomicRanges::GRanges(win_ranges$V1,
                               ranges=IRanges::IRanges(win_ranges$start,
                                                       win_ranges$end))
  return(gr)
}
get_genomic_range <- function(grs, cds, win) {
  end1 <- as.numeric(as.character(GenomicRanges::end(grs[win])))
  end2 <- as.numeric(as.character(GenomicRanges::start(grs[win])))
  win_range <- cds[(fData(cds)$bp1 < end1 &
                                fData(cds)$bp1 > end2) |
                               (fData(cds)$bp2 < end1 &
                                  fData(cds)$bp2 > end2), ]
  win_range <-
    win_range[as.character(fData(win_range)$chr) ==
                gsub("chr", "",
                     as.character(GenomicRanges::seqnames(grs[win]))),]
  fData(win_range)$mean_bp <-
    (as.numeric(as.character(fData(win_range)$bp1)) +
       as.numeric(as.character(fData(win_range)$bp2)))/2

  return(win_range)
}





#run test


get_bagged_matrices=function(seuratobj){

  exc<-seuratobj
  DefaultAssay(exc)<-"ATAC"
  exc_atac_mat<-GetAssayData(exc)


  exc_cds<-create_atac_cds(exc_atac_mat)
  exc_coords<-get_atac_coords(exc_cds)
  # exc_bag_cds<-make_cicero_cds(exc_cds,exc_coords)

  DefaultAssay(exc)<-"RNA"
  exc_rna_mat<-GetAssayData(exc)


  exc_atac_cds<-create_atac_cds(exc_atac_mat)
  exc_bag_atac_cds<-create_bagged_pseuduobulk(exc_atac_cds,tsne_dims=exc_coords)

  exc_rna_cds<-create_exp_cds(exc_rna_mat)
  exc_bag_rna_cds<-create_bagged_pseuduobulk(exc_rna_cds,tsne_dims=exc_coords)



  final_list<-list(exc_bag_atac_cds,exc_bag_rna_cds)
  names(final_list)<-c("atac","rna")
  return(final_list)



}

merged_seurat<-readRDS("/mnt/data/ALEX/PROJECTS/MULTIOME_PD/CICERO/04OCT_ANNOTATED_combined_multiome_seuratobj.rds")
merged_seurat<-SplitObject(merged_seurat,split.by="RNAannotation")
exc<-merged_seurat$Excitatory
DefaultAssay(exc)<-"ATAC"
exc_atac_mat<-GetAssayData(exc)


exc_cds<-create_atac_cds(exc_atac_mat)
exc_coords<-get_atac_coords(exc_cds)
exc_bag_cds<-make_cicero_cds(exc_cds,exc_coords)

DefaultAssay(exc)<-"RNA"
exc_rna_mat<-GetAssayData(exc)


exc_atac_cds<-create_atac_cds(exc_atac_mat)
exc_bag_atac_cds<-create_bagged_pseuduobulk(exc_atac_cds,tsne_dims=exc_coords)

exc_rna_cds<-create_exp_cds(exc_rna_mat)
exc_bag_rna_cds<-create_bagged_pseuduobulk(exc_rna_cds,tsne_dims=exc_coords)



rna_mat<-log2(cpm(exc_bag_rna_cds)+1)
atac_mat<-log2(cpm(exc_bag_atac_cds)+1)



cor.test(atac_mat["chr7_23205804_23206687",],rna_mat["GPNMB",],method="spearman")



correlation_wrap=function(seuratobj,ld_window,gene,)
