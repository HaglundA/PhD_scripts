




source("/rds/general/user/ah3918/home/R_FUNCTIONS/MRC_coloc_funcs.r")

outdir="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/HIP/"

library(GenomicFeatures)
library(GenomicRanges)
library(MatrixEQTL)
library(TreeQTL)
library(Seurat)
library(SNPRelate)
library(SeqArray)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(org.Hs.eg.db)
library(dplyr)

run_id="standardized"
cellQTL(inputfile=paste0(outdir,"HIP_seuratObject_Aug.rds"),
               Step1=TRUE,
               project_name="HIP",
               method="agg_cpm",
               Step2=TRUE,
               gds_inputfile="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/meqtl_files/MRC_60samples_imputed_SNV_genotyping_Sep.gds",
               outputdir=outdir,
               filter_pseudobulk_thresh=3,
               covadj_matrixeqtl=TRUE,
               run_id=run_id,
               cov_file="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/meqtl_files/MRC_covariates_allsamples.txt",
               covs=c("Age","PMI","Gender","Sample_Source"),
               filter_genotype_matrix=TRUE,pval_thresh_cis=1,geno_filter_type = "2.2")



mateqtldir=paste0(outdir,"MatrixEQTL_IO/",Sys.Date(),"_",run_id,"/")
mateqtlouts<-eqtl_rds(mateqtldir=mateqtldir)
message("Saving mateqtlouts..")
saveRDS(mateqtlouts,paste0(mateqtldir,Sys.Date(),"_full_mateqtlouts.rds"))
