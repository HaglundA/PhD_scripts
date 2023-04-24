




chunk_chroms=function(snppositions,windowsize){
    
    chrs=unique(snppositions$chrom)
    
    
    allchunks=data.frame()
    
    for(i in 1:length(chrs)){
        chr=chrs[i]
        message(paste0("Chunking ",chr, "..."))
        snpstokeep=geno_loc[geno_loc$chrom %in% chrs[i],]
        start=snpstokeep$position[1]
        end=start+windowsize

        window=data.frame(start=start,end=end,chr=chr)
        nchunks=1
        while(end<snpstokeep$position[nrow(snpstokeep)]){

            start=start+windowsize
            end=end+windowsize

            window=rbind(window,c(start,end,chr))
            nchunks=nchunks+1

        }
        

        message(paste0("Chromosome split into ",nchunks," chunks."))
        allchunks=rbind(allchunks,window)
        allchunks$start=as.numeric(allchunks$start)
        allchunks$end=as.numeric(allchunks$end)

    }
    return(allchunks)
}




run_atlasqtl_wrap=function(geno_mat,
exp_list,
windowsize=2e6,
geno_loc,
ref_chunks=NULL,
specify_chr=NULL,
exp_loc){

    #first chunk geno_loc

    if(!is.null(ref_chunks)){
        chunks=ref_chunks
    }else{
     chunks=chunk_chroms(geno_loc,windowsize=windowsize)
    }

    if(!is.null(specify_chr)){
        chunks=chunks[chunks$chr %in% specify_chr,]
    }

   
    reslist=list()
    for(i in 1:nrow(chunks)){
        start_time=Sys.time()

        chr=chunks$chr[i]
        start=chunks$start[i]-1
        end=chunks$end[i]+1
        snpstokeep=geno_loc %>% 
        mutate(position=as.numeric(position)) %>%
        filter(chrom %in% chr) %>%
        filter(position>start & position<end)

        chrom=chunks$chr[i]
        ##filter genes
        genes=filter(exp_loc,chr %in% chrom)

        genes=genes %>%
        filter((left>start & left<end) | (right>start &left<end))


        snpstokeep=snpstokeep$annot
        genestokeep=genes$geneid

        message(paste0(length(snpstokeep)," SNPs and ",length(genestokeep)," genes used in the analysis."))

        input_test=lapply(exp_list,function(x){
    
        geno_mat=geno_mat[complete.cases(geno_mat),]
        common_names=intersect(colnames(x),colnames(geno_mat))
        
        snps=geno_mat[snpstokeep,common_names]
    #     snps=mutate_all(snps,function(x) as.numeric(x))
                
        snps=t(snps)
        expr=x[genestokeep,common_names]
        expr=t(expr)
        
        return(list(expr=expr,snps=snps))
        })

        res=run_atlasqtl(input_test)
        reslist[[i]]=res
        end_time=Sys.time()
        print(difftime(end_time,start_time))
        


    }
    return(reslist)

}


run_atlasqtl=function(inputlist){


    RNGkind("L'Ecuyer-CMRG")
    seed <- 1
    set.seed(seed)

    n_cpus <- 8

    mu_t <- 1 # prior number of SNPs expected to be associated with each expression level
    sd_t <- 2 # prior standard deviation for this number
    start_time<-Sys.time()

    rt=system.time(list_out_test <- parallel::mclapply(names(inputlist), function(type) {

    snps <- inputlist[[type]]$snps
    expr <- inputlist[[type]]$expr

    stopifnot(rownames(snps) == rownames(expr))
    
    atlasqtl(Y = expr, X = snps, 
            p0 = c(mu_t, sd_t^2), 
            add_collinear_back = TRUE)
    
    }, mc.cores = n_cpus))

    return(list_out_test)
}