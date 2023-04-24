


tmp<-vector()
for(i in 2:8){
  y<-1*10^-i
  x<-5*10^-i
  tmp<-c(tmp,x,y)
}

basedir="/mnt/data/ALEX/COLOC_GENOMEWIDE_RUN/MatrixEQTL_IO/HYPRCOLOC/"

tmp1<-vector()
for(i in 1:length(tmp)){
  outdir=paste0(basedir,"GWAS_Threshold_",as.character(tmp[i]))
  setwd(outdir)
  finaldf<-readRDS("finaldf_hits.rds")
  hits<-finaldf[which(finaldf$posterior_prob>0),]

  tmp1<-c(tmp1,nrow(hits))


}
tmp2<-data.frame(hits=tmp1,thresh=tmp)
write.table(tmp2,paste0(basedir,"All_hits_all_thresholds.txt"))
