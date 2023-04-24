



append_hg38_snps=function(hg19GWAS,outfile="conversiontable.csv",MAGMA=T){



  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
  library(org.Hs.eg.db)

  message("Reading in GWAS..")
  hg19GWAS<-read.table(hg19GWAS,header=TRUE,stringsAsFactors=F)

  message("Creating GRanges object..")
  chromlocs<-as.character(hg19GWAS$SNP)
  chromlocs<-strsplit(chromlocs,":")
  chromlocs<-do.call(rbind.data.frame,chromlocs)
  colnames(chromlocs)<-c("Chr","Pos")

  chronly<-strsplit(as.character(chromlocs$Chr),"r")
  chronly<-do.call(rbind.data.frame,chronly)

  final<-data.frame(chrom=chronly[,2],position=chromlocs$Pos)
  # final<-final[order(final$chrom),]
  final$paste<-paste0(final$chrom,":",final$position)
  hg19GWAS_gr<-makeGRangesFromDataFrame(final,
  seqnames.field="chrom",
  start.field="position",
  end.field="position")

  seqlevelsStyle(hg19GWAS_gr)="UCSC"

  message("Performing liftover..")

  chain<-import.chain("/mnt/data/ALEX/Tools/hg19ToHg38.over.chain")
  # chain<-import.chain("/rds/general/user/ah3918/home/Tools/hg19ToHg38.over.chain")

  hg38GWAS_gr<-liftOver(hg19GWAS_gr,chain)
  hg38GWAS_gr<-unlist(hg38GWAS_gr)

  seqlevelsStyle(hg19GWAS_gr)="NCBI"
  seqlevelsStyle(hg38GWAS_gr)="NCBI"

  hg19snp<-SNPlocs.Hsapiens.dbSNP144.GRCh37
  hg19rsids<-snpsByOverlaps(hg19snp,hg19GWAS_gr)

  hg38snp<-SNPlocs.Hsapiens.dbSNP151.GRCh38
  hg38rsids<-snpsByOverlaps(hg38snp,hg38GWAS_gr)

  hg19_rsids_df<-data.frame(rsid=as.data.frame(hg19rsids)$RefSNP_id,
                         chr_pos=paste0("chr",as.data.frame(hg19rsids)$seqnames,":",as.data.frame(hg19rsids)$pos))

  hg38_rsids_df<-data.frame(rsid=as.data.frame(hg38rsids)$RefSNP_id,
                         chr_pos=paste0("chr",as.data.frame(hg38rsids)$seqnames,":",as.data.frame(hg38rsids)$pos))

  hg19GWAS_df<-as.data.frame(hg19GWAS_gr)
  hg38GWAS_df<-as.data.frame(hg38GWAS_gr)

  hg19GWAS_df<-hg19GWAS_df[rownames(hg19GWAS_df) %in% rownames(hg38GWAS_df),]

  conversiontable<-data.frame(hg19_pos=paste0("chr",hg19GWAS_df$seqnames,":",hg19GWAS_df$start),
                                       hg38_pos=paste0("chr",hg38GWAS_df$seqnames,":",hg38GWAS_df$start))


  conversiontable$hg38_chr_only<-hg38GWAS_df$seqnames
  conversiontable$hg38_pos_only<-hg38GWAS_df$start

  test<-hg38_rsids_df[match(conversiontable$hg38_pos,hg38_rsids_df$chr_pos),]
  conversiontable$rsid_hg38<-test$rsid
  conversiontable$hg19_pos<-as.character(conversiontable$hg19_pos)
  conversiontable$hg38_pos<-as.character(conversiontable$hg38_pos)
  conversiontable$rsid_hg38<-as.character(conversiontable$rsid_hg38)



  n_snps_before<-nrow(hg19GWAS)

  hg19GWAS<-hg19GWAS[match(conversiontable$hg19_pos,hg19GWAS$SNP),]
  hg19GWAS$rsid_hg38<-conversiontable$rsid_hg38
  hg19GWAS$hg38_pos<-conversiontable$hg38_pos
  hg19GWAS$hg38_chr_only<-conversiontable$hg38_chr_only
  hg19GWAS$hg38_pos_only<-conversiontable$hg38_pos_only

  hg19GWAS<-hg19GWAS[complete.cases(hg19GWAS),]

  n_snps_after<-nrow(hg19GWAS)

  snps_lost<-n_snps_before - n_snps_after

  message(paste0("Conversion complete. ",snps_lost," SNPs lost over conversion"))

  if(MAGMA==T){

    magmadf<-data.frame(SNP=hg19GWAS$rsid_hg38,
      CHR=hg19GWAS$hg38_chr_only,
      BP=hg19GWAS$hg38_pos_only,
      P=hg19GWAS$p,
      NOBS=hg19GWAS$N_cases)

    write.table(magmadf,"06APR21.txt",row.names=F,quote=F,sep="\t")
  } else {
    write.table(hg19GWAS, outfile)
  }


}



liftover_hg19_hg38=function(GWAS){
  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
  library(org.Hs.eg.db)

  message("Reading in GWAS..")
  gwas_name<-GWAS


  GWAS<-read.table(GWAS,header=TRUE,stringsAsFactors=F)

  if(length(grep("beta",colnames(GWAS)))==0 && length(grep("OR",colnames(GWAS)))==1){
    GWAS$beta<-log(GWAS$OR)
  }

  if(length(grep("chr",GWAS$chr))==0){
    GWAS$chr<-paste0("chr",GWAS$chr)
  }

  chromlocs<-cbind(GWAS$chr,GWAS$pos)
  colnames(chromlocs)<-c("Chr","Pos")
  chromlocs<-as.data.frame(chromlocs)





  # message("Creating GRanges object..")
  # chromlocs<-as.character(GWAS$SNP)
  # chromlocs<-strsplit(chromlocs,":")
  # chromlocs<-do.call(rbind.data.frame,chromlocs)

    chronly<-strsplit(as.character(chromlocs$Chr),"r")
    chronly<-do.call(rbind.data.frame,chronly)

   final<-data.frame(chrom=chronly[,2],position=chromlocs$Pos)
   final<-final[complete.cases(final),]

  final$paste<-paste0(final$chrom,":",final$position)
  GWAS_gr<-makeGRangesFromDataFrame(final,
  seqnames.field="chrom",
  start.field="position",
  end.field="position")

  seqlevelsStyle(GWAS_gr)="UCSC"

  message("Performing liftover..")

  chain<-import.chain("/mnt/data/ALEX/Tools/hg19ToHg38.over.chain")
  # chain<-import.chain("/rds/general/user/ah3918/home/Tools/hg19ToHg38.over.chain")

  hg38GWAS_gr<-liftOver(GWAS_gr,chain)
  hg38GWAS_gr<-unlist(hg38GWAS_gr)

  seqlevelsStyle(GWAS_gr)="NCBI"
  seqlevelsStyle(hg38GWAS_gr)="NCBI"

  hg19snp<-SNPlocs.Hsapiens.dbSNP144.GRCh37
  hg19rsids<-snpsByOverlaps(hg19snp,GWAS_gr)

  hg38snp<-SNPlocs.Hsapiens.dbSNP151.GRCh38
  hg38rsids<-snpsByOverlaps(hg38snp,hg38GWAS_gr)
  #
  hg19_rsids_df<-data.frame(rsid=as.data.frame(hg19rsids)$RefSNP_id,
                         chr_pos=paste0("chr",as.data.frame(hg19rsids)$seqnames,":",as.data.frame(hg19rsids)$pos))
  #
  hg38_rsids_df<-data.frame(rsid=as.data.frame(hg38rsids)$RefSNP_id,
                         chr_pos=paste0("chr",as.data.frame(hg38rsids)$seqnames,":",as.data.frame(hg38rsids)$pos))
  #
  GWAS_df<-as.data.frame(GWAS_gr)
  hg38GWAS_df<-as.data.frame(hg38GWAS_gr)

  GWAS_df<-GWAS_df[rownames(GWAS_df) %in% rownames(hg38GWAS_df),]
  #
  conversiontable<-data.frame(hg19_pos=paste0("chr",GWAS_df$seqnames,":",GWAS_df$start),
                                       hg38_pos=paste0("chr",hg38GWAS_df$seqnames,":",hg38GWAS_df$start))


  conversiontable$hg38_chr_only<-hg38GWAS_df$seqnames
  conversiontable$hg38_pos_only<-hg38GWAS_df$start

  test<-hg38_rsids_df[match(conversiontable$hg38_pos,hg38_rsids_df$chr_pos),]
  conversiontable$rsid_hg38<-test$rsid
  conversiontable$hg19_pos<-as.character(conversiontable$hg19_pos)
  conversiontable$hg38_pos<-as.character(conversiontable$hg38_pos)
  conversiontable$rsid_hg38<-as.character(conversiontable$rsid_hg38)

  n_snps_before<-nrow(GWAS)
  message(paste0("Started with ",n_snps_before," variants."))

  GWAS<-GWAS[match(conversiontable$rsid_hg38,GWAS$rsid),]
  GWAS$rsid_hg38<-conversiontable$rsid_hg38
  GWAS$hg38_pos<-conversiontable$hg38_pos
  GWAS$hg38_chr_only<-conversiontable$hg38_chr_only
  GWAS$hg38_pos_only<-conversiontable$hg38_pos_only

  GWAS<-GWAS[complete.cases(GWAS),]

  n_snps_after<-nrow(GWAS)

  snps_lost<-n_snps_before - n_snps_after

  GWAS$pos<-GWAS$hg38_pos_only
  message(paste0(snps_lost," variants were lost from conversion."))

  write.table(GWAS,paste0("liftedover/",gwas_name,"_liftedover.txt"))
}
