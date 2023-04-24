


# CONVERT GENE IDs

# library(gdsfmt)
# library(SNPRelate)
# library(SeqArray)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

convert_geneids=function(genelist,
  format=c("entrezID","ENSEMBL","SYMBOL"),
  conversion=c("ENSEMBL","EntrezID","SYMBOL")){


  library(org.Hs.eg.db)

  suppressMessages(symbols <- mapIds(org.Hs.eg.db, keys = genelist, keytype = format, column=conversion))

  symbols<-as.data.frame(symbols)
  symbols$pre_conversion<-rownames(symbols)
  genesbefore<-nrow(symbols)
  # symbols<-symbols[complete.cases(symbols),]

  genesafter<-genesbefore-nrow(symbols)

  if(genesafter>0){
    message(paste0(genesafter," genes were lost/didn't match Ensembl IDs. Double-check gene names!"))
    return(symbols)
  } else {
    return(symbols)
  }


}
