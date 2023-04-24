






run_MR=function(pruned_region,trait){

  gwas_beta<-as.data.frame(pruned_region$betas)$GWAS
  gwas_se<-as.data.frame(pruned_region$std.errors)$GWAS

  betas<-as.data.frame(pruned_region$betas)[,trait]
  ses<-as.data.frame(pruned_region$std.errors)[,trait]
  MRInputObject<-mr_input(bx=betas,bxse=ses,by=gwas_beta,byse=gwas_se)


}


library(MendelianRandomization)

biglist_se<-readRDS("biglist_se.rds")
biglist_beta<-readRDS("biglist_beta.rds")

gwas_se<-as.data.frame(biglist_se[[25]])$GWAS
gwas_beta<-as.data.frame(biglist_beta[[25]])$GWAS

betas<-as.data.frame(biglist_beta[[25]])
betas<-

ses<-as.data.frame(biglist_se[[25]])
ses<-




 MRInputObject<-mr_input(bx=betas,bxse=ses,by=gwas_beta,byse=gwas_se)

 ivwobj<-mr_ivw(MRInputObject,
    model = "default",
    robust = FALSE,
    penalized = FALSE,
    correl = FALSE,
    weights = "simple",
    psi = 0,
    distribution = "normal",
    alpha = 0.05)


medianobj<-mr_median(MRInputObject,
  weighting = "weighted",
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265)

EggerObject<-mr_egger(MRInputObject,
robust = FALSE,
penalized = FALSE,
correl = FALSE,
distribution = "normal",
alpha = 0.05)

MaxLikObject <- mr_maxlik(MRInputObject,
model = "default",
correl = FALSE,
psi = 0,
distribution = "normal",
alpha = 0.05)




harmonize_effects=function(mateqtlouts,geno_maf,gwas){

  #make sure GWAS has alleles


}
