


manhattan_plot=function(coloc_data_path,
                       trait,
                       ldlinkr_token=NULL,
                       snp_loc_path,
                       pop="CEU",
                        leadsnp=NULL,
                       build="grch38"){
    
    
    
    message("Reading in files..")
    coloc_data=readRDS(coloc_data_path)
    snp_locs=as.data.frame(data.table::fread(snp_loc_path))


    cols=lapply(coloc_data[[1]],function(x){
      return(colnames(x))
    })

    region_index=grep(trait,cols)

    #extract pvalues
    pvaldf=coloc_data$pvalues[[region_index]]

    #The LD matrix calculation accepts 1000 SNPs max, so we'll take the top 1000 SNPs based on p-value (we are not interested in correlation between weak eQTLs)
    pvaldf=pvaldf[order(pvaldf[,trait]),]
    pvaldf=pvaldf[1:1000,]

    snps=rownames(pvaldf)
    snps=snps[snps %in% snp_locs$annot]

    #get LD matrix. It only accepts 1000 SNPs max, so we'll take the top 1000 SNPs based on p-value (we are not interested in correlation between weak eQTLs)

    ldmat=LDlinkR::LDmatrix(snps,
                            pop=pop,
                            genome_build=build,
                            token=ldlinkr_token)

    rownames(ldmat)=ldmat$RS_number
    ldmat$RS_number=NULL

    ##not all SNPs get used in the LD mat. Filter.
    pvaldf=pvaldf[rownames(ldmat),]
    pval_trait=pvaldf[,trait]
    pval_gwas=pvaldf[,"GWAS"]
    snps=rownames(pvaldf)

    markersdf=data.frame(marker=snps,
                       chr=snp_locs[match(snps,snp_locs$annot),]$chrom,
                       pos=snp_locs[match(snps,snp_locs$annot),]$position)

    markersdf$chr=as.data.frame(do.call(rbind,strsplit(markersdf$chr,"r")))$V2

    zdf=data.frame(Trait1=qnorm(pval_trait),Trait2=qnorm(pval_gwas))
    colnames(zdf)=c(trait,"GWAS")
    
    if(!is.null(leadsnp)){
        
        gassocplot2::stack_assoc_plot(markersdf,zdf,ldmat,traits=colnames(zdf),top.marker = leadsnp)
    }else{
        gassocplot2::stack_assoc_plot(markersdf,zdf,ldmat,traits=colnames(zdf))
        
    }
 
    
}



##example usage 

manhattan_plot(coloc_data_path="../AD_2022_controlonly/coloc_data.rds",
               trait="Microglia.BIN1",
               snp_loc_path="../../MatrixEQTL_IO/snp_chromlocations.csv",
             ldlinkr_token="f3d054e6c0ee")




manhattan_plot_simple=function(coloc_data_path,
                       trait,
                       snp_loc_path){


    message("Reading in files..")
    coloc_data=readRDS(coloc_data_path)
    snp_locs=as.data.frame(data.table::fread(snp_loc_path))


    cols=lapply(coloc_data[[1]],function(x){
      return(colnames(x))
    })

    region_index=grep(trait,cols)

    #extract pvalues
    pvaldf=coloc_data$pvalues[[region_index]]

    pval_trait=pvaldf[,trait]
    pval_gwas=pvaldf[,"GWAS"]
    snps=rownames(pvaldf)

    markersdf=data.frame(marker=snps,
                       chr=snp_locs[match(snps,snp_locs$annot),]$chrom,
                       pos=snp_locs[match(snps,snp_locs$annot),]$position)

    markersdf$chr=as.data.frame(do.call(rbind,strsplit(markersdf$chr,"r")))$V2

    zdf=data.frame(Trait1=pval_trait,Trait2=pval_gwas)
    colnames(zdf)=c(trait,"GWAS")

    markersdf=cbind(markersdf,zdf)
    manhattan.plot<-function(chr, pos, pvalue,
        sig.level=NA, annotate=NULL, ann.default=list(),
        should.thin=T, thin.pos.places=2, thin.logp.places=2,
        xlab="Chromosome", ylab=expression(-log[10](p-value)),
        col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {

        if (length(chr)==0) stop("chromosome vector is empty")
        if (length(pos)==0) stop("position vector is empty")
        if (length(pvalue)==0) stop("pvalue vector is empty")

        #make sure we have an ordered factor
        if(!is.ordered(chr)) {
            chr <- ordered(chr)
        } else {
            chr <- chr[,drop=T]
        }

        #make sure positions are in kbp
        if (any(pos>1e6)) pos<-pos/1e6;

        #calculate absolute genomic position
        #from relative chromosomal positions
        posmin <- tapply(pos,chr, min);
        posmax <- tapply(pos,chr, max);
        posshift <- head(c(0,cumsum(posmax)),-1);
        names(posshift) <- levels(chr)
        genpos <- pos + posshift[chr];
        getGenPos<-function(cchr, cpos) {
            p<-posshift[as.character(cchr)]+cpos
            return(p)
        }

        #parse annotations
        grp <- NULL
        ann.settings <- list()
        label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5,
            col=NULL, fontface=NULL, fontsize=NULL, show=F)
        parse.label<-function(rawval, groupname) {
            r<-list(text=groupname)
            if(is.logical(rawval)) {
                if(!rawval) {r$show <- F}
            } else if (is.character(rawval) || is.expression(rawval)) {
                if(nchar(rawval)>=1) {
                    r$text <- rawval
                }
            } else if (is.list(rawval)) {
                r <- modifyList(r, rawval)
            }
            return(r)
        }

        if(!is.null(annotate)) {
            if (is.list(annotate)) {
                grp <- annotate[[1]]
            } else {
                grp <- annotate
            }
            if (!is.factor(grp)) {
                grp <- factor(grp)
            }
        } else {
            grp <- factor(rep(1, times=length(pvalue)))
        }

        ann.settings<-vector("list", length(levels(grp)))
        ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

        if (length(ann.settings)>1) {
            lcols<-trellis.par.get("superpose.symbol")$col
            lfills<-trellis.par.get("superpose.symbol")$fill
            for(i in 2:length(levels(grp))) {
                ann.settings[[i]]<-list(pch=pch,
                    col=lcols[(i-2) %% length(lcols) +1 ],
                    fill=lfills[(i-2) %% length(lfills) +1 ],
                    cex=cex, label=label.default);
                ann.settings[[i]]$label$show <- T
            }
            names(ann.settings)<-levels(grp)
        }
        for(i in 1:length(ann.settings)) {
            if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
            ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label,
                parse.label(ann.settings[[i]]$label, levels(grp)[i]))
        }
        if(is.list(annotate) && length(annotate)>1) {
            user.cols <- 2:length(annotate)
            ann.cols <- c()
            if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
                ann.cols<-match(names(annotate)[-1], names(ann.settings))
            } else {
                ann.cols<-user.cols-1
            }
            for(i in seq_along(user.cols)) {
                if(!is.null(annotate[[user.cols[i]]]$label)) {
                    annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label,
                        levels(grp)[ann.cols[i]])
                }
                ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]],
                    annotate[[user.cols[i]]])
            }
        }
        rm(annotate)

        #reduce number of points plotted
        if(should.thin) {
            thinned <- unique(data.frame(
                logp=round(-log10(pvalue),thin.logp.places),
                pos=round(genpos,thin.pos.places),
                chr=chr,
                grp=grp)
            )
            logp <- thinned$logp
            genpos <- thinned$pos
            chr <- thinned$chr
            grp <- thinned$grp
            rm(thinned)
        } else {
            logp <- -log10(pvalue)
        }
        rm(pos, pvalue)
        gc()

        #custom axis to print chromosome names
        axis.chr <- function(side,...) {
            if(side=="bottom") {
                panel.axis(side=side, outside=T,
                    at=((posmax+posmin)/2+posshift),
                    labels=levels(chr),
                    ticks=F, rot=0,
                    check.overlap=F
                )
            } else if (side=="top" || side=="right") {
                panel.axis(side=side, draw.labels=F, ticks=F);
            }
            else {
                axis.default(side=side,...);
            }
        }

        #make sure the y-lim covers the range (plus a bit more to look nice)
        prepanel.chr<-function(x,y,...) {
            A<-list();
            maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
            A$ylim=c(0,maxy);
            A;
        }

        xyplot(logp~genpos, chr=chr, groups=grp,
            axis=axis.chr, ann.settings=ann.settings,
            prepanel=prepanel.chr, scales=list(axs="i"),
            panel=function(x, y, ..., getgenpos) {
                if(!is.na(sig.level)) {
                    #add significance line (if requested)
                    panel.abline(h=-log10(sig.level), lty=2);
                }
                panel.superpose(x, y, ..., getgenpos=getgenpos);
                if(!is.null(panel.extra)) {
                    panel.extra(x,y, getgenpos, ...)
                }
            },
            panel.groups = function(x,y,..., subscripts, group.number) {
                A<-list(...)
                #allow for different annotation settings
                gs <- ann.settings[[group.number]]
                A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]
                A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
                A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
                A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
                A$x <- x
                A$y <- y
                do.call("panel.xyplot", A)
                #draw labels (if requested)
                if(gs$label$show) {
                    gt<-gs$label
                    names(gt)[which(names(gt)=="text")]<-"labels"
                    gt$show<-NULL
                    if(is.character(gt$x) | is.character(gt$y)) {
                        peak = which.max(y)
                        center = mean(range(x))
                        if (is.character(gt$x)) {
                            if(gt$x=="peak") {gt$x<-x[peak]}
                            if(gt$x=="center") {gt$x<-center}
                        }
                        if (is.character(gt$y)) {
                            if(gt$y=="peak") {gt$y<-y[peak]}
                        }
                    }
                    if(is.list(gt$x)) {
                        gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
                    }
                    do.call("panel.text", gt)
                }
            },
            xlab=xlab, ylab=ylab,
            panel.extra=panel.extra, getgenpos=getGenPos, ...
        );
    }

    test=list(
        manhattan.plot(chr=markersdf$chr,pos=markersdf$pos,pvalue=markersdf$GWAS,xlab="GWAS"),
        manhattan.plot(chr=markersdf$chr,pos=markersdf$pos,pvalue=markersdf[,trait],xlab=trait))
    return(grid.arrange(test[[1]],test[[2]],nrow=2))

}

    