preprocess_data <- function(qdnaseq_data, pileup_data, phased_bcf, sample_map, sex, build='hg19', max_phaseable_distance=20000, min_total_count_for_baf=5, min_allele_count_for_het=10) {
    #max_phaseable_distance <- 20000; min_total_count_for_baf <- 5; min_allele_count_for_het <- 10

    all_chrs <- c(1:22,'X','Y','MT')
    if(sex=='female') {
        diploid_chrs <- c(1:22,'X')
    } else {
        diploid_chrs <- c(1:22)       
    }

    ## load 1Mb bin counts from QDNAseq. We need to annotate these bins with allele-specific counts (for bins containing het-SNPs)
    ## this probably means the het-SNPs bins need to perfectly align with the QDNAseq bins
    cnt <- readRDS(qdnaseq_data)
    map <- fread(sample_map)
    cnt <- merge(cnt, map[,c('qdnaseq_name','proper_name'),with=F], by.x='file', by.y='qdnaseq_name', all.x=T)
    cnt <- cnt[!is.na(proper_name)]
    cnt[,file:=NULL]
    setnames(cnt,c('feature','chromosome','start','end','count','count_raw','binsize_kbp','proper_name'),c('bin','chr','bin_start','bin_end','count','count_raw','binsize_kbp','sample'))
    cnt$chr <- factor(cnt$chr, levels=all_chrs)
    cnt <- cnt[!is.na(chr),]
    cnt <- cnt[order(chr,bin_start,bin_end),]
    cnt[,bin_length:=(bin_end-bin_start+1)/1e6]
    cnt$count_raw <- as.numeric(cnt$count_raw)
    samples <- map$proper_name

    ## For chrM, correct bin counts for variable bin length
    ## nb: bin_length was normalized to have max,median of 1
    #cnt[chr=='MT', count_raw:=count_raw/bin_length] 
    #cnt[chr=='MT', count:=count/bin_length]

    ## load phased SNPs
    message('Loading phased SNP bcf-file: ',phased_bcf,' ...')
    bcf <- scanBcf(phased_bcf)
    phased_snps <- data.table(CHROM=bcf$CHROM, POS=bcf$POS, ID=bcf$ID, REF=bcf$REF, ALT=bcf$ALT, PHASE=bcf$GENO$GT[,1])

    ## subset for chrs that are diploid in normal given sex
    phased_snps$CHROM <- gsub('chr','',phased_snps$CHROM)
    phased_snps <- phased_snps[CHROM %in% diploid_chrs,]
    phased_snps[,POS2:=POS+1]
    phased_snps$CHROM <- factor(phased_snps$CHROM, levels=diploid_chrs) ## use all chrs for levels
    setkey(phased_snps,'CHROM','POS','POS2')

    ## annotate the phased snps with 1Mb-binned counts
    bins <- cnt[!duplicated(bin),c('bin','chr','bin_start','bin_end'),with=F]
    setkey(bins,'chr','bin_start','bin_end')

    phased_snps <- foverlaps(phased_snps, bins, type='within')
    phased_snps <- phased_snps[!is.na(bin),]
    phased_snps$bin <- factor(phased_snps$bin, levels=unique(phased_snps$bin))
    phased_snps$bin.i <- as.integer(phased_snps$bin)

    ## define phasing blocks
    message('Defining phasing blocks with maximum phaseable distance: ',max_phaseable_distance,' ...')
    phased_snps <- block_phased_snps(phased_snps, max_phaseable_distance=max_phaseable_distance)

    ## get the number of SNPs per phasing block in each bin. Choose the largest block in each bin and discard the rest.
    collapse_blocks <- function(phased_snps) {
        block_start = min(phased_snps$POS)
        block_end = max(phased_snps$POS)
        block_snps = length(unique(phased_snps$POS))
        list(block_start=block_start, block_end=block_end, block_snps=block_snps)
    }
    blocks <- phased_snps[, collapse_blocks(.SD), by=c('CHROM','bin','block')]
    tmp <- blocks[order(bin, block_snps, decreasing=T),]
    tmp <- tmp[!duplicated(bin),]
    blocks <- tmp[order(block),]
    blocks[,block_midpoint:=round((block_start+block_end)/2)]

    ## annotate each SNP position with %GC
    message('Annotating each SNP position with %GC ...')
    get_GCbias <- function(i, phased_snps) { 
        gc <- data.table(CHROM=i, POS=phased_snps[CHROM==i,(POS)], GC=getGCpct(i, phased_snps[CHROM==i,(POS)], gbuild=build))
    }
    gcbias <- rbindlist(lapply(1:22, get_GCbias, phased_snps))
    phased_snps$GC <- gcbias$GC
    phased_snps$GC <- round(phased_snps$GC, 2)

    ## load allele-specific read counts for the phased het-SNPs for each sample
    message('Loading snp-pileup data: ',pileup_data,' ...')
    d <- fread(pileup_data, sep=',')
    d$Chromosome <- gsub('chr','',d$Chromosome)
    d <- d[Chromosome %in% diploid_chrs]
    d$Chromosome <- factor(d$Chromosome, levels=diploid_chrs)
    d <- merge(d, phased_snps[,c('CHROM','POS','PHASE','bin','block','GC'),with=F], by.x=c('Chromosome','Position'), by.y=c('CHROM','POS'), all.x=T)
    d <- d[!is.na(block),]
    front_fields <- names(d)[grepl('File',names(d))==F]
    back_fields <- names(d)[!names(d) %in% front_fields]
    d <- d[,c(front_fields, back_fields),with=F]

    ## correct raw read count for each allele and sample for GC%
    tmp <- data.table(pos=1:ncol(d), orig_name=names(d))
    tmp[grepl('File', orig_name),pileup_pos:=as.integer(gsub('File','',strtrim(orig_name,nchar(orig_name)-1)))]
    tmp[grepl('File', orig_name),allele:=substr(orig_name,nchar(orig_name), nchar(orig_name))]
    tmp <- merge(tmp, map[,c('pileup_pos','proper_name'),with=F], by='pileup_pos', all.x=T)
    tmp <- tmp[order(pos),]
    tmp[!is.na(pileup_pos), new_name:=paste0(proper_name,'.',allele)] 
    tmp[is.na(pileup_pos), new_name:=orig_name]
    names(d) <- tmp$new_name 
    
    message('Correcting read counts for each SNP allele for GC% bias (done for each sample individually)')
    field_names <- tmp[allele %in% c('R','A') & !is.na(pileup_pos),(new_name)]
    for(f in field_names) {
        message(f)
        x <- d[,c('GC',f),with=F]
        names(x)[2] <- 'count'
        qs <- quantile(x$count,c(0.005,0.995),na.rm=T)
        x[(count < qs[1] | count > qs[2]), count:=NA]
        grand_mean <- mean(x$count, na.rm=T)
        x[,diff_from_mean:=count - grand_mean]
        fit <- lm(diff_from_mean ~ GC, data=x)
        predicted_diff_from_mean <- predict(fit, newdata=d, type='response')
        d[[f]] <- d[[f]]-predicted_diff_from_mean
    }

    ## extract the corrected read counts per block for each allele
    count10_A <- d[PHASE=='1|0',c('block',paste0(samples,'.A')),with=F]; names(count10_A) <- c('block',samples)
    count01_R <- d[PHASE=='0|1',c('block',paste0(samples,'.R')),with=F]; names(count01_R) <- c('block',samples)
    count1 <- rbind(count10_A, count01_R)
    count01_A <- d[PHASE=='0|1',c('block',paste0(samples,'.A')),with=F]; names(count01_A) <- c('block',samples)
    count10_R <- d[PHASE=='1|0',c('block',paste0(samples,'.R')),with=F]; names(count10_R) <- c('block',samples)
    count2 <- rbind(count01_A, count10_R)
    collapse_block <- function(count) {
        out <- colSums(count)     
        as.list(out)
    }
    block_counts1 <- count1[,collapse_block(.SD),by=c('block')]
    block_counts2 <- count2[,collapse_block(.SD),by=c('block')]

    ## generate the BAF matrix
    message('Generating BAF matrix for blocked SNPs in each 1Mb bin ...')
    block_counts <- merge(block_counts1, block_counts2, by='block', all=F) ## merge the corrected counts for each allele per block. Only retain blocks with counts for both alleles
    baf <- matrix(nrow=nrow(block_counts), ncol=length(samples))
    rownames(baf) <- block_counts$block
    colnames(baf) <- samples
    for(s in samples) {
        baf[,s] <- block_counts[[paste0(s,'.x')]] / (block_counts[[paste0(s,'.x')]] + block_counts[[paste0(s,'.y')]])
    }
    bad_blocks <- which(is.na(rowSums(baf)))
    if(length(bad_blocks) > 0) baf <- baf[-bad_blocks, ]
    baf <- cbind(block=as.integer(rownames(baf)), as.data.table(baf))

    ## Multiply the binned-count by the largest block's BAF to get allele-specific counts for this bin.
    bafM <- data.table::melt(baf, id.var='block')
    names(bafM) <- c('block','sample','baf')
    bafM$sample <- as.character(bafM$sample)
    bafM <- merge(bafM, blocks[,c('block','bin'),with=F], by='block', all.x=T)
    bafM <- bafM[!is.na(bin),]
    d <- merge(cnt, bafM, by=c('sample','bin'), all.x=T)
    d <- d[order(sample,chr, bin_start, bin_end, block),]
    d[is.nan(baf),baf:=NA]
    d[!is.na(baf),count1:=count * baf]
    d[!is.na(baf),count2:=count * (1-baf)]
    d[is.na(baf), count1:=count]
    d[is.na(baf), count2:=0]
    d[is.na(baf),Ref:='.']
    d[is.na(baf),Alt:='.']
    d[!is.na(baf),Ref:='N']
    d[!is.na(baf),Alt:='N']
    d[,position:=round((bin_start + bin_end)/2)]
    #browser()
    d <- d[,c('chr','position','Ref','Alt','sample','count1','count2','count_raw'),with=F]
    d <- d[order(chr,position,sample),]

    message('Adding logR and BAF values for each bin ...')
    d[is.na(count_raw),count_raw:=0]
    d[,dp:=count1+count2]
    d[Ref=='N', BAF:=count1 / dp]
    d[BAF < 0 | BAF > 1, BAF:=NA]

    d[,bin:=paste0(chr,':',position)]
    d$bin <- as.integer(factor(d$bin, levels=unique(d$bin)))
    get_logR <- function(d) {
        mid <- median(d$dp[d$chr %in% c(1:22,'X')],na.rm=T)
        d$logR <- log2(d$dp / mid)
        d
    }
    d <- d[,get_logR(.SD),by=sample]

    message('Adding chr-arms ...')
    arms <- get_chr_arms()$arms
    arms[,arm_start:=arm_start*1e6]  
    arms[,arm_end:=arm_end*1e6]  
    setkey(arms,'chr','arm_start','arm_end')
    d[,start:=(position/1e6) - 0.5]
    d[,end:=(position/1e6) + 0.5]
    setkey(d,'chr','start','end')
    d <- foverlaps(d, arms, type='within')
    setnames(d,c('chr','position'),c('Chromosome','Position'))
    d[Chromosome=='MT', arm:='']
    d <- d[!is.na(arm),]
    d[,c('start','end','arm_start','arm_end'):=NULL]
    d[,charm:=paste0(Chromosome,arm)]
    message('Done!')
    d
}





block_phased_snps <- function(phased_snps, max_phaseable_distance) {
    ## within each 1Mb bin, we create phasing blocks 
    require(Rcpp)

    src <-
        "NumericVector rcpp_phaseblock(NumericVector pos, NumericVector bin, double max_phaseable_distance){

        // initialize variables    
        int n = pos.length();
        NumericVector block (n); // create a vector for each row's block value

        // get starting values before loop
        double currentbin = bin[0]; // starting bin
        double currentpos = pos[0]; // starting position
        double blockcount = 0;
        block[0] = 0;

        // loop over all positions for this bin
        for(int i=1; i<n; ++i){
            double nextpos = pos[i]; // new pos
            double nextbin = bin[i]; // new bin
            double posdiff = nextpos - currentpos; // distance from the previous pos

            if(posdiff > max_phaseable_distance || nextbin > currentbin) {
                blockcount = blockcount + 1; // increment the block count
            }

            // note down the block number for this position
            block[i] = blockcount; 

            // update the pos and bin variables
            currentpos = nextpos;
            currentbin = nextbin;
        }
        return(block);
    }"

    cppFunction(src)
    pos <- phased_snps$POS
    bin.i <- phased_snps$bin.i
    blocks <- rcpp_phaseblock(pos, bin.i, max_phaseable_distance)
    phased_snps$block <- blocks
    phased_snps
}

write_distance_matrix <- function (dm, file) {
    write.table(dm, file = file, sep = "\t", quote = FALSE, col.names = NA)
}

get_logR <- function(pu, pl, na, nb) log2((2*(1-pu) + pu*(na + nb)) / (2*(1-pu) + pu*pl))

get_BAF <- function(pu, pl, na, nb) (1 - pu + pu*nb) / (2 - 2*pu + pu*(na+nb))

write_tsv <- function (d, file, sep = "\t", quote = F, row.names = F, ...) {
    write.table(d, file = file, sep = sep, quote = quote, row.names = row.names, ...)
}

d2m <- function (dt) {
    rows <- dt[[1]]
    if ("data.table" %in% class(dt)) {
        dt <- dt[, c(2:ncol(dt)), with = F]
    }
    else if (class(dt) == "data.frame") {
        dt <- dt[, c(2:ncol(dt))]
    }
    else {
        stop("enter a data.table or a data.frame")
    }
    m <- as.matrix(dt)
    rownames(m) <- rows
    m
}


get_na_nb <- function(p, t, x) {
    ## I think this works?
    r <- x$logR
    b <- x$BAF
    g <- x$germline_copies
    c <- 2^r
    c1 <- 2^(r+1)

    # this is slightly incorrect because the median normal ploidy can have GC other than 2 (median should always be 2)
    #na = (b*g*p*c - b*g*c - b*p*c*t - g*p*c + g*p + g*c - g + p*c*t - p + 1)/p
    #nb = -1*(b*g*p*c - b*g*c - b*p*c*t - p + 1) / p

    ## NB: here we model the median admixed ploidy as 2*(1-p)+p*t
    ## this works because the logR values are calculated as log2(x/median(x))
    na = -1*( b*p*c*t - b*p*c1 + b*c1 - g*p + g - p*c*t + p*c1 + p - c1 - 1 ) / p
    nb = (b*p*c*t - b*p*c1 + b*c1 + p - 1) / p        

    ## go over previous values where BAF was NA
    NAs <- is.na(b)
    tmp_na <- (g*(p-1) + c*(p*(t-2)+2)) / p
    na[NAs] <- tmp_na[NAs]
    nb[NAs] <- as.numeric(NA)

    x$na <- na
    x$nb <- nb
    x
}



get_loglik <- function(pu, pl, obj) {

    sample_dat <- obj$sample_dat
    sample_seg <- obj$sample_seg
    logR_CDF <- obj$logR_CDF
    BAF_CDF <- obj$BAF_CDF

    sample_seg2 <- get_na_nb(p=pu, t=pl, sample_seg)
    sample_seg2 <- sample_seg2[!is.na(na) & !is.na(nb) & !is.na(BAF) & !is.na(logR),]
    sample_seg2$na <- round(sample_seg2$na)
    sample_seg2$nb <- round(sample_seg2$nb)
    sample_seg2[na < 0, na:=0]
    sample_seg2[nb < 0, nb:=0]

    sample_seg2$logR_expected <- get_logR(pu, pl, sample_seg2$na, sample_seg2$nb)
    sample_seg2$BAF_expected <- get_BAF(pu, pl, sample_seg2$na, sample_seg2$nb)
    sample_dat2 <- merge(sample_dat, sample_seg2[,c('segment','logR_expected','BAF_expected'),with=F], by='segment', all.x=T)
    sample_dat2 <- sample_dat2[!is.na(logR_expected) & !is.na(BAF_expected) & Chromosome %in% c(1:22,'X'),] ## do not use Y or MT for loglik 

    ## get the expected logR and BAF values if na,nb are the nearest integers
    diff_from_logR_expected <- sample_dat2$i.logR - sample_dat2$logR_expected
    diff_from_logR_expected <- diff_from_logR_expected[!is.na(diff_from_logR_expected)]
    diff_from_BAF_expected <- sample_dat2$i.BAF - sample_dat2$BAF_expected
    diff_from_BAF_expected <- diff_from_BAF_expected[!is.na(diff_from_BAF_expected)]

    ## calculate the summed loglikelihood of the bins' observed logRs and BAFs, given the nearest integer copy numbers
    loglik_logR <- log(c(1-logR_CDF(diff_from_logR_expected[diff_from_logR_expected >= 0]), logR_CDF(diff_from_logR_expected[diff_from_logR_expected < 0])))
    min_valid_logR_loglik <- min(loglik_logR[!is.infinite(loglik_logR)])
    max_valid_logR_loglik <- max(loglik_logR[!is.infinite(loglik_logR)])
    loglik_logR[loglik_logR==-Inf] <- min_valid_logR_loglik
    loglik_logR[loglik_logR==Inf] <- max_valid_logR_loglik

    loglik_BAF <- log(c(1-BAF_CDF(diff_from_BAF_expected[diff_from_BAF_expected >= 0]), BAF_CDF(diff_from_BAF_expected[diff_from_BAF_expected < 0])))
    min_valid_BAF_loglik <- min(loglik_BAF[!is.infinite(loglik_BAF)])
    max_valid_BAF_loglik <- max(loglik_BAF[!is.infinite(loglik_BAF)])
    loglik_BAF[loglik_BAF==-Inf] <- min_valid_BAF_loglik
    loglik_BAF[loglik_BAF==Inf] <- max_valid_BAF_loglik

    loglik <- sum(loglik_logR) + sum(loglik_BAF)
    loglik
}




get_fit <- function(obj, dipLogR=NA, purity=NA, ploidy=NA, max_homdel=1e8, max_neg=0, min_ploidy=1, max_ploidy=8, cores=1, best_only=T) {
    require(parallel)

    sample_name <- obj$samplename
    sample_dat <- obj$sample_dat
    sample_seg <- obj$sample_seg
    logR_CDF <- obj$logR_CDF
    loglik_BAF <- obj$loglik_BAF
    fit_chrs <- obj$fit_chrs
    message('Fitting for chrs: ', paste(fit_chrs, collapse=', '))
    sample_dat <- sample_dat[Chromosome %in% fit_chrs]
    sample_seg <- sample_seg[Chromosome %in% fit_chrs]

    if(!is.na(dipLogR) & (is.na(purity) | is.na(ploidy))) {        
        message('dipLogR provided')
        auto <- T
        purity <- seq(0.05,1,by=0.001)
        ploidy <- (2*(purity + 2^(-dipLogR) - 1)) / purity
        myfits <- data.table(pu=purity, pl=round(ploidy,3))
        myfits$dipLogR <- dipLogR
    } else if(!is.na(purity) & !is.na(ploidy)) {
        auto <- F
        message('Purity/ploidy provided')
        myfits <- data.table(pu=purity, pl=ploidy)
        myfits$dipLogR = round(log2(( 2*(1-purity) + purity*2 ) / ( 2*(1-purity) + purity*ploidy )), 4)
    } else {
        auto <- T
        message('grid-searching purity/ploidy combinations')
        purity <- seq(0.05,1,by=0.025)
        #ploidy <- seq(1,6,by=0.025)
        ploidy <- seq(min_ploidy, max_ploidy, by=0.025)
        myfits <- as.data.table(expand.grid(pu=purity, pl=ploidy))
        myfits$dipLogR = round(log2(( 2*(1-myfits$pu) + myfits$pu*2 ) / ( 2*(1-myfits$pu) + myfits$pu*myfits$pl )), 2)
    }

    ## for all fits, first filter out impossible fits based on negative copies and too much homozygous deletion
    qc_fit <- function(i, myfits, sample_seg) {
        this.pu <- myfits$pu[i]
        this.pl <- myfits$pl[i]
        this.dipLogR <- myfits$dipLogR[i]
        fit <- get_na_nb(this.pu, this.pl, sample_seg)
        len_neg <- sum(sample_seg$seg_length[fit$na <= -0.5 | fit$nb <= -0.5],na.rm=T)
        len_homdel <- sum(sample_seg$seg_length[round(fit$na) == 0 & round(fit$nb) == 0],na.rm=T)
        len_na <- sum(sample_seg$seg_length[is.na(fit$na) | is.na(fit$nb)])
        list(pu=this.pu, pl=this.pl, dipLogR=this.dipLogR, len_neg=len_neg, len_homdel=len_homdel, len_na=len_na)
    }
    n_fits <- nrow(myfits)
    profile_list <- mclapply(1:n_fits, qc_fit, myfits, sample_seg, mc.cores=cores) 
    myfits <- rbindlist(profile_list)
    if(auto==T) myfits <- myfits[len_neg <= max_neg & len_homdel <= max_homdel & pl <= max_ploidy,]

    get_possible_fits_for_sample <- function(myfits, obj, cores) {
        try_fits <- function(i, myfits, obj) {
            pu <- myfits$pu[i]
            pl <- myfits$pl[i]
            loglik <- get_loglik(pu, pl, obj)
            list(pl=pl, pu=pu, loglik=loglik)
        }
        n_fits <- nrow(myfits)
        l <- mclapply(1:n_fits, try_fits, myfits, obj, mc.cores=cores)
        res <- rbindlist(l)
        myfits$loglik <- res$loglik
        myfits
    }

    if(nrow(myfits) > 0) {
        myfits <- get_possible_fits_for_sample(myfits, obj, cores)
        myfits <- myfits[order(loglik,decreasing=T),]
        if(best_only==T) myfits <- myfits[1,]
        myfits$sample <- sample_name
        myfits    
    } else {
        NULL
    }
}




plot_fit <- function(fit, obj, int_copies=F, highlight_seg=c()) {
    require(cowplot)
    require(ggplot2)
    sample_dat <- obj$sample_dat
    sample_seg <- obj$sample_seg
    sample_seg[,global_seg_start:=(seg_start/1e6) + global_start_mb]
    sample_seg[,global_seg_end:=(seg_end/1e6) + global_start_mb]
    sample_seg[, highlight:=F]
    sample_seg[segment %in% highlight_seg, highlight:=T]
    sname <- paste0(obj$samplename,' (purity=',fit$pu,', ploidy=',fit$pl,', loglik=',round(fit$loglik,3),', dipLogR=',fit$dipLogR,')')
    sex <- obj$sex
    highlight <- sample_seg[highlight==T,]
    valid_chrs <- unique(sample_seg[germline_copies>0,(Chromosome)])
    chr <- get_chr_arms()$chr

    plot_sample_seg <- sample_seg[Chromosome %in% valid_chrs]
    plot_sample_dat <- sample_dat[Chromosome %in% valid_chrs]
    plot_chr <- chr[chr %in% valid_chrs]
    plot_chr$chr <- factor(plot_chr$chr, levels=valid_chrs)
    plot_sample_dat$chr <- factor(plot_sample_dat$chr, levels=valid_chrs)
    plot_sample_seg$chr <- factor(plot_sample_seg$chr, levels=valid_chrs)
    centerline <- median(plot_sample_seg$logR,na.rm=T)

    p1 <- ggplot(plot_sample_dat, aes(x=global_pos, y=i.logR)) + 
        geom_point(color='#bfbfbf',size=1.25,pch=16) +
        geom_hline(yintercept=centerline, color='green', linewidth=0.25) +
        geom_hline(yintercept=fit$dipLogR, color='purple', linewidth=0.25) +
        geom_vline(xintercept=c(0,plot_chr$global_end)) +
        #scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2,2,by=1)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr, expand=c(0,0)) +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start,xend=global_seg_end,y=logR,yend=logR),color='blue',linewidth=0.75,lineend='round') +
        theme_fit(base_size=12) +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x=NULL, y='logR', subtitle=sname)
    p2 <- ggplot(plot_sample_dat, aes(x=global_pos, y=i.BAF)) + 
        geom_point(color='#bfbfbf',size=1.25,pch=16) +
        geom_vline(xintercept=c(0,plot_chr$global_end)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr,expand=c(0,0)) +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start,xend=global_seg_end,y=BAF,yend=BAF),color='blue',linewidth=0.75,lineend='round') +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start,xend=global_seg_end,y=1-BAF,yend=1-BAF),color='blue',linewidth=0.75,lineend='round') +
        theme_fit(base_size=12) +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x=NULL, y='BAF')

    ascn <- get_na_nb(p=fit$pu, t=fit$pl, x=plot_sample_seg)
    ascn[,tcn:=na+nb]
    ascn[na >= nb,mcn:=nb]
    ascn[na < nb,mcn:=na]
    ascn[is.na(nb) & !is.na(na), tcn:=na]
    ascn[is.na(nb) & !is.na(na), mcn:=NA]
    ascn[tcn >= 10, tcn:=log10(tcn)+9]
    ascn[mcn >= 10, tcn:=log10(mcn)+9]
    ascn[,tcn_ad_from_int:=abs(tcn-round(tcn))] 
    ascn[,mcn_ad_from_int:=abs(mcn-round(mcn))] 
    ascn[tcn < -1, tcn:=-1]
    ascn[mcn < -1, mcn:=-1]

    if(int_copies==T) {
        ascn[,tcn:=round(tcn)]
        ascn[,mcn:=round(mcn)]
    }

    p3 <- ggplot(ascn) + 
        scale_y_continuous(breaks=c(seq(0,10,by=2),log10(100)+9), labels=c(seq(0,10,by=2),100)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr, expand=c(0,0)) 
    if(nrow(highlight) > 0) p3 <- p3 + annotate("rect", xmin=highlight$global_seg_start, xmax=highlight$global_seg_end, ymin=-0.5, ymax=log10(100)+9, fill='steelblue', alpha=0.1)
    p3 <- p3 + 
        geom_hline(yintercept=seq(0,10,by=1), color='#bfbfbf', linewidth=0.25) +
        geom_vline(xintercept=c(0,plot_chr$global_end),linewidth=0.25) +
        geom_segment(aes(x=global_seg_start,xend=global_seg_end,y=tcn,yend=tcn,alpha=tcn_ad_from_int),linewidth=1,color='black',lineend='round') +
        geom_segment(aes(x=global_seg_start,xend=global_seg_end,y=mcn,yend=mcn,alpha=mcn_ad_from_int),linewidth=1,color='red',lineend='round') +
        scale_alpha(range=c(0.5,0), limits=c(0,0.75), breaks=seq(0,0.5,by=0.25), name='AD from integer CN') +
        theme_fit(base_size=12) +
        theme(legend.position='bottom') +
        labs(x='Genomic posititon', y='Copy number') 
    p <- plot_grid(p1, p2, p3, align='v', ncol=1, rel_heights=c(1.05,1,1.2), axis='lr')
    p
}


get_obj_for_sample <- function(this.sample, seg, d, sex) { 
    sample_seg <- seg[sample==this.sample]
    
    ## add the expected number of germline copies given sex
    #sample_seg[Chromosome %in% 1:22, germline_copies:=2]
    #sample_seg[Chromosome %in% 'MT', germline_copies:=NA]
    #if(sex=='female') {
    #    sample_seg[Chromosome %in% 'X', germline_copies:=2]
    #    sample_seg[Chromosome %in% 'Y', germline_copies:=0]
    #} else {
    #    sample_seg[Chromosome %in% 'X', germline_copies:=1]
    #    sample_seg[Chromosome %in% 'Y', germline_copies:=1]   
    #}
    fit_chrs <- unique(sample_seg[germline_copies>0,(Chromosome)])

    sample_dat <- d[sample==this.sample]
    sample_dat[,Position2:=Position+1]
    setkey(sample_dat,'Chromosome','Position','Position2')
    toadd <- sample_seg[,c('Chromosome','seg_start','seg_end','BAF','logR','segment','global_start_mb','germline_copies'),with=F]
    setkey(toadd,'Chromosome','seg_start','seg_end')
    sample_dat <- foverlaps(sample_dat, toadd, type='within')
    sample_dat$Chromosome <- factor(sample_dat$Chromosome, levels=c(1:22,'X','Y','MT'))
    sample_dat <- sample_dat[order(Chromosome, Position),]
    sample_dat[,global_pos:=Position/1e6+global_start_mb]
    sample_dat[,c('seg_start','seg_end'):=NULL]

    ## in each segment, what does the distribution around the logR and BAF estimates look like?
    qc <- sample_dat[Chromosome %in% fit_chrs, c('segment','i.BAF','i.logR','BAF','logR'),with=F]
    qc[,diff_from_BAF:=i.BAF - BAF]
    qc[,diff_from_logR:=i.logR - logR]
    logR_CDF <- ecdf(qc$diff_from_logR)
    BAF_CDF <- ecdf(qc$diff_from_BAF)
    list(sample_dat=sample_dat, sample_seg=sample_seg, logR_CDF=logR_CDF, BAF_CDF=BAF_CDF, sex=sex, samplename=this.sample, fit_chrs=fit_chrs)
}



theme_fit <- function (base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_bw(base_size = base_size, base_line_size = base_line_size,
        base_rect_size = base_rect_size) %+replace% theme(line = element_line(colour = "black", linewidth = base_line_size, linetype = 1, lineend = "round"),
        text = element_text(colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = F), 
        axis.text = element_text(colour = "black", size = rel(0.8)), 
        axis.ticks = element_line(colour = "black", size = rel(1)), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(), #line(colour = "black", size = rel(1)), 
        legend.key = element_blank()) 
} 


run_multipcf <- function(d, normal_sample, sex, build, penalty=70, tmpdir='.', cleanup=T, seed=as.integer(Sys.time())) {
    if(!is.na(seed)) set.seed(42)

    if(sex=='female') {
        valid_chr <- c(1:22,'X')
        gender <- 'XX'
    } else {
        gender <- 'XY'
        valid_chr <- c(1:22,'X','Y')
    }
    d <- d[Chromosome %in% valid_chr,]
    d$Chromosome <- factor(d$Chromosome, levels=valid_chr)
    
    tumor_samples <- unique(d[sample!=normal_sample,(sample)])

    message('Preparing to run multipcf from ASCAT.')
    t_logR <- data.table::dcast(Chromosome + Position ~ sample, data=d[sample %in% tumor_samples], value.var='logR')
    t_logR <- t_logR[order(Chromosome, Position),]
    n_logR <- dcast(Chromosome + Position ~ sample, data=d[sample %in% normal_sample,], value.var='logR')
    n_logR <- n_logR[order(Chromosome, Position),]
    t_BAF <- data.table::dcast(Chromosome + Position ~ sample, data=d[sample %in% tumor_samples], value.var='BAF')
    t_BAF <- t_BAF[order(Chromosome, Position),]
    n_BAF <- dcast(Chromosome + Position ~ sample, data=d[sample %in% normal_sample,], value.var='BAF')
    n_BAF <- n_BAF[order(Chromosome, Position),]

    if(!dir.exists(tmpdir)) dir.create(tmpdir)
    tumor_logR_file <- file.path(tmpdir,'Tumor_LogR.txt')
    tumor_BAF_file <- file.path(tmpdir,'Tumor_BAF.txt')
    germline_logR_file <- file.path(tmpdir,'Germline_LogR.txt')
    germline_BAF_file <- file.path(tmpdir,'Germline_BAF.txt')
    write.table(t_logR, file = tumor_logR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_logR, file = germline_logR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(t_BAF, file = tumor_BAF_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_BAF, file = germline_BAF_file, sep = "\t", quote = FALSE, col.names = NA)

    ascat.bc = ascat.loadData(Tumor_LogR_file = tumor_logR_file,
                              Tumor_BAF_file = tumor_BAF_file,
                              Germline_LogR_file = germline_logR_file,
                              Germline_BAF_file = germline_BAF_file,
                              gender = rep(gender,nrow(t_logR), genomeVersion = build)) 
    message('Running ascat.asmultipcf() with penalty=',penalty,', seed=',seed,'.')
    ascat.bc = ascat.asmultipcf(ascat.bc, out.dir=NA, seed=seed, penalty=penalty, refine=F)

    if(cleanup==T) {
        message('Cleaning up after multipcf.')
        trash <- file.remove(c(tumor_logR_file,germline_logR_file,tumor_BAF_file,germline_BAF_file))
    }
    
    ## extract segmented BAF and logR data from ASCAT
    get_extracted <- function(field, ascat.bc) {
        tmp <- ascat.bc[[field]]
        extract_sample_data <- function(i, tmp) {
            x <- tmp[[i]]
            this.sample <- colnames(x)
            out <- data.table(bin=rownames(x), value=as.numeric(x))
            out$bin <- factor(out$bin, levels=out$bin)
            out$sample <- this.sample
            out
        }
        l <- lapply(1:length(tmp), extract_sample_data, tmp)    
        bafs <- rbindlist(l)
        bafs
    }
    baf <- get_extracted('Tumor_BAF_segmented', ascat.bc)
    names(baf) <- c('bin','BAF','sample')
    logR <- ascat.bc[['Tumor_LogR_segmented']]
    logR <- cbind(bin=rownames(logR), as.data.table(logR))
    logR <- data.table::melt(logR, id.var='bin')
    names(logR) <- c('bin','sample','logR')
    snps <- ascat.bc$SNPpos
    snps <- cbind(bin=rownames(snps), as.data.table(snps))
    multipcf.extracted <- merge(baf, logR, by=c('sample','bin'), all=T)
    multipcf.extracted <- merge(multipcf.extracted, snps, by='bin', all.multipcf.extracted=T)
    multipcf.extracted$Chromosome <- factor(multipcf.extracted$Chromosome, levels=valid_chr)
    multipcf.extracted <- multipcf.extracted[order(sample, Chromosome, Position),]
    multipcf.extracted$bin <- as.integer(gsub('bin','',multipcf.extracted$bin)) 
    multipcf.extracted$Chromosome <- factor(multipcf.extracted$Chromosome, levels=unique(multipcf.extracted$Chromosome))
    multipcf.extracted <- multipcf.extracted[order(Chromosome, Position, sample),]
    multipcf.extracted
}


get_chr_arms <- function(chr_lengths_file=here('~/lab_repos/lpASCN/misc/chr_lengths_b37.txt'), cytoBand_file=here('~/lab_repos/lpASCN/misc/cytoBand_b37.txt')) {

    ## get a table of hg19 chromosome lengths with global start/end/midpoint positions (in Mb) 
    chr <- fread(chr_lengths_file)
    chr$length <- chr$length / 1e6
    chr$chr_start <- 0
    chr$chr_end <- chr$length
    chr$global_start <- rep(as.integer(NA),nrow(chr))
    chr$global_end <- rep(as.integer(NA),nrow(chr))
    chr$global_midpoint <- rep(as.integer(NA),nrow(chr))
    chr$global_start[1] <- 0
    chr$global_end[1] <- chr$length[1]
    chr$global_midpoint[1] <- (chr$chr_start[1] + chr$chr_end[1])/2
    for(i in 2:nrow(chr)) {
        chr$global_start[i] <- sum(chr$length[1:(i-1)])
        chr$global_end[i] <- chr$global_start[i]+chr$length[i]
        chr$global_midpoint[i] <- (chr$global_start[i] + chr$global_end[i])/2
    }
    chr$chr <- factor(chr$chr, levels=c(1:22,'X','Y','MT'))
    setkey(chr,'chr','chr_start','chr_end')

    bands <- fread(cytoBand_file)
    names(bands) <- c('chr','band_start','band_end','region','stain')
    bands[,arm:=strtrim(region,1)]
    bands[,chr:=gsub('chr','',chr)]
    bands[,c('stain','region'):=NULL]
    collapse_arm <- function(bands) {
        arm_start <- min(bands$band_start) / 1e6
        arm_end <- max(bands$band_end) / 1e6
        list(arm_start=arm_start, arm_end=arm_end) 
    }
    arms <- bands[,collapse_arm(.SD),by=c('chr','arm')]
    arms$chr <- factor(arms$chr, levels=c(1:22,'X','Y','MT'))    
    arms <- merge(arms, chr, by='chr', all=T)
    arms[arm=='p', arm_start:=chr_start] 
    arms[arm=='q', arm_end:=chr_end] 
    #arms[,global_arm_start:=global_start + arm_start]
    #arms[,global_arm_end:=global_start + arm_end]
    arms <- arms[,c('chr','arm','arm_start','arm_end'),with=F]
    arms[chr=='MT', arm:='']
    arms[chr=='MT',arm_start:=0.0]
    arms[chr=='MT',arm_end:=16569/1e6]
    arms <- arms[!is.na(chr),]
    setkey(arms,'chr','arm_start','arm_end')
    list(chr=chr, arms=arms)
}


get_segments <- function(d, normal_sample, sex, build='hg19', penalty=70, seed=42, tmpdir='.', cleanup=T, binsize=1e6) {
    #build <- 'hg19'; penalty <- 70; cleanup <- T; seed=42; normal_sample <- 'Normal1'; sex <- 'female'; tmpdir <- '.'
   
    all_chrs <- c(1:22,'X','Y','MT')
    if(sex=='female') {
        diploid_chrs <- c(1:22,'X')
        haploid_chrs <- c()
    } else {
        diploid_chrs <- c(1:22)
        haploid_chrs <- c('X','Y')       
    }
    l <- run_multipcf(d, normal_sample=normal_sample, sex=sex, build=build, penalty=penalty, cleanup=cleanup, seed=seed, tmpdir=tmpdir)

    ## add chrM back in
    chrM <- d[Chromosome=='MT',]
    toadd <- chrM[sample!=normal_sample,c('bin','sample','BAF','logR','Chromosome','Position'),with=F]
    l <- rbind(l, toadd)

    ## annotate with chromosome arms 
    chr_arms <- get_chr_arms()
    arms <- chr_arms$arms
    arms[,chr_arm:=paste0(chr,arm)]
    l[,pos_mb1:=Position/1e6]
    l[,pos_mb2:=(Position+1)/1e6]
    setkey(l,'Chromosome','pos_mb1','pos_mb2')
    setkey(arms,'chr','arm_start','arm_end')
    l <- foverlaps(l, arms, type='within') 
    l <- l[!is.na(bin) & !is.na(arm),]

    ## extract logR and BAF matrices
    l$Chromosome <- factor(l$Chromosome, levels=unique(l$Chromosome))
    logR <- data.table::dcast(bin + Chromosome + Position + chr_arm + arm + arm_start + arm_end ~ sample, value.var='logR', data=l)
    bins <- logR[,c('bin','Chromosome','Position','chr_arm','arm','arm_start','arm_end'),with=F]
    bin_arm <- logR[,c('bin','arm'),with=F]
    logR[,c('arm','Chromosome','Position','chr_arm','arm_start','arm_end'):=NULL]
    logR <- d2m(logR)
    BAF <- d2m(data.table::dcast(bin ~ sample, value.var='BAF', data=l))

    ## get aligned segments across all samples based on the segmented BAF and logR values and the chromosome arm
    n <- nrow(BAF)
    this.logR <- logR[1,]
    this.BAF <- BAF[1,]
    this.arm <- bin_arm$arm[1]
    segment <- rep(as.integer(NA), n)
    segment[1] <- 1
    segment.count <- 1
    for(i in 2:n) {
        next.logR <- logR[i,]
        next.BAF <- BAF[i,]
        next.arm <- bin_arm$arm[i]
        valid.logR <- !is.na(next.logR) & !is.na(this.logR)
        valid.BAF <- !is.na(next.BAF) & !is.na(this.BAF)
        if( any(next.logR[valid.logR]!=this.logR[valid.logR]) | any(next.BAF[valid.BAF]!=this.BAF[valid.BAF]) | next.arm!=this.arm ) {
            segment.count <- segment.count + 1
        }
        this.logR <- next.logR; this.BAF <- next.BAF; this.arm <- next.arm;
        segment[i] <- segment.count + 1 
    }
    bins$segment <- segment
    #browser()
    ll <- merge(l, bins[,c('bin','segment'),with=F], by='bin', all.x=T)  
    
    ## get wide segment-level data for logR and BAF
    collapse_segment <- function(ll) {
        seg_start <- min(ll$Position)
        seg_end <- max(ll$Position)
        BAF <- unique(ll$BAF[!is.na(ll$BAF)])
        logR <- unique(ll$logR[!is.na(ll$logR)])
        n_bins <- length(unique(ll$bin))
        if(length(BAF) > 1 | length(BAF) > 1) stop('Multiple values per segment?')
        list(seg_start=seg_start, seg_end=seg_end, BAF=BAF, logR=logR, n_bins=n_bins)
    }
    seg <- ll[,collapse_segment(.SD), by=c('sample','Chromosome','arm','segment','arm_start','arm_end')]
    chr <- chr_arms$chr
    seg <- merge(seg, chr[,c('chr','global_start'),with=F], by.x='Chromosome', by.y='chr', all.x=T)
    setnames(seg,'global_start','global_start_mb')
    seg[,arm_start:=arm_start*1e6]
    seg[,arm_end:=arm_end*1e6]
    seg[,seg_start:=seg_start-(binsize/2)]
    seg[,seg_end:=seg_end+(binsize/2)]
    seg[seg_start < arm_start, seg_start:=arm_start] 
    seg[seg_end > arm_end, seg_end:=arm_end] 
    seg[,seg_length:=(seg_end-seg_start)]
    
    norm <- seg[!duplicated(segment),]
    norm[,logR:=NULL]
    norm[,BAF:=NULL]
    norm$Chromosome <- factor(norm$Chromosome, levels=c(1:22,'X','Y','MT'))
    extract_normal_for_segment <- function(this.segment, norm, d) {
        message(this.segment)
        start <- norm[segment==this.segment,(seg_start)] 
        end <- norm[segment==this.segment,(seg_end)] 
        this.chr <- norm[segment==this.segment,(Chromosome)] 
        logR <- median(d[sample==normal_sample & Position >= start & Position <= end & Chromosome==this.chr, (logR)], na.rm=T)
        BAF <- median(d[sample==normal_sample & Position >= start & Position <= end & Chromosome==this.chr, (BAF)], na.rm=T)
        list(segment=this.segment, logR=logR, BAF=BAF)
    } 
    segments <- unique(norm$segment)
    norm_dat <- rbindlist(lapply(segments, extract_normal_for_segment, norm, d))
    norm <- merge(norm, norm_dat, by='segment', all.x=T)
    norm$sample <- normal_sample
    
    out <- rbind(norm, seg)
    out[Chromosome %in% 1:22, germline_copies:=2]
    out[Chromosome=='MT', germline_copies:=NA]
    if(sex=='female') {
        out[Chromosome %in% 'X', germline_copies:=2]
        out[Chromosome %in% 'Y', germline_copies:=0]
    } else {
        out[Chromosome %in% 'X', germline_copies:=1]
        out[Chromosome %in% 'Y', germline_copies:=1]
    }
    out
}


refine_segments <- function(d, seg, k.max, min.bins.for.clustering=10, normal_sample, sex, build, penalty, seed=42, tmpdir='.', cleanup=T, binsize=1e6) {
    require(Ckmeans.1d.dp)

    ## combine the preprocessed data with a first-round of segments
    d[,Position2:=Position+1]
    seg <- seg[!duplicated(segment),c('segment','Chromosome','seg_start','seg_end','n_bins'),with=F]
    setkey(d,'Chromosome','Position','Position2')
    setkey(seg,'Chromosome','seg_start','seg_end')
    d <- foverlaps(d, seg, type='any')

    cluster_baf <- function(dd, k.max) {
        dd <- dd[order(BAF,decreasing=F),]
        clus <- Ckmeans.1d.dp(dd$BAF, k=c(1,k.max))
        dd$k <- clus$cluster
        k.max.used <- max(clus$cluster)
        if(k.max.used > 2)  dd[!k %in% c(1,k.max.used), BAF:=NA]
        dd
    }

    valid_clusters <- seg[n_bins >= min.bins.for.clustering, (segment)]
    d1 <- d[is.na(segment) | !segment %in% valid_clusters]
    d1$k <- as.integer(NA)
    d2 <- d[segment %in% valid_clusters, cluster_baf(.SD,k.max), by=c('sample','segment')]
    dnew <- rbind(d1, d2)
    dnew$Chromosome <- factor(dnew$Chromosome, levels=c(1:22,'X','Y','MT'))
    dnew <- dnew[order(Chromosome, Position, sample),]

    ## now regenerate the segments
    segnew <- get_segments(dnew, normal_sample=normal_sample, sex=sex, build=build, penalty=penalty, seed=seed, tmpdir='.', cleanup=cleanup, binsize=binsize)
    list(d=dnew, seg=segnew)
}


get_aligned_segs <- function(mb, standard, seg) {
    ## subset for unique solutions
    same_values <- grep('same_',names(mb),value=T)
    mb$sameID <- apply(mb[,(same_values),with=F], 1, paste, collapse='')
    mb <- mb[!duplicated(sameID),]
    m_values <- grep('m_',names(mb),value=T)
    m_values <- data.table::melt(mb[,c(m_values,'w'),with=F], id.var='w')
    m_values$variable <- as.character(m_values$variable)
    m_values$variable <- gsub('m_','',m_values$variable)
    setnames(m_values,'variable','sample')
    setnames(m_values,'value','m')

    b_values <- grep('b_',names(mb),value=T)
    b_values <- data.table::melt(mb[,c(b_values,'w'),with=F], id.var='w')
    b_values$variable <- as.character(b_values$variable)
    b_values$variable <- gsub('b_','',b_values$variable)
    setnames(b_values,'variable','sample')
    setnames(b_values,'value','b')

    mb_values <- merge(m_values, b_values, by=c('w','sample'))
    toadd <- data.table(w=(seq(1:nrow(mb))-1), sample=standard, m=1, b=0)
    mb_values <- rbind(mb_values, toadd)

    same_values <- grep('same_',names(mb),value=T)
    same_values <- data.table::melt(mb[,c(same_values,'w'),with=F], id.var='w')
    same_values$variable <- as.character(same_values$variable)
    same_values$variable <- as.integer(gsub('same_','',same_values$variable))
    setnames(same_values,'variable','segment')
    setnames(same_values,'value','same')
    same_values <- d2m(data.table::dcast(segment ~ w, value.var='same', data=same_values))
    solutions <- apply(same_values, 2, paste, collapse='')
    toadd <- data.table(segment=as.integer(rownames(same_values)), same_values[,!duplicated(solutions)])
    names(toadd) <- c('segment', paste0('sol',1:(ncol(toadd)-1)))

    seg <- merge(seg, mb_values[w==0,], by=c('sample'), all.x=T)
    seg <- seg[!is.na(w),]
    seg <- seg[order(sample, segment),]
    seg[,logRadj:=logR*m - b]
    seg[,global_seg_start:=seg_start/1e6 + global_start_mb]
    seg[,global_seg_end:=seg_end/1e6 + global_start_mb]
    seg <- merge(seg, toadd, by='segment', all.x=T)
    seg$standard <- standard
    seg
} 


get_alignment_plot <- function(aligned_seg, solution) {
    standard <- unique(aligned_seg$standard)
    sol <- aligned_seg[[solution]]
    same_values <- unique(aligned_seg[!is.na(sol) & sol==1, (segment)])

    tmp1 <- aligned_seg[Chromosome %in% c(1:22,'X','Y'),c('sample','segment','global_seg_start','global_seg_end','logR'),with=F]
    tmp1$type <- 'Original'
    tmp2 <- aligned_seg[Chromosome %in% c(1:22,'X','Y'),c('sample','segment','global_seg_start','global_seg_end','logRadj'),with=F]
    setnames(tmp2,'logRadj','logR')
    tmp2$type <- 'Aligned'
    plotdat <- rbind(tmp1, tmp2)
    plotdat$same <- plotdat$segment %in% same_values
    plotdat$type <- factor(plotdat$type, levels=c('Original','Aligned'))
    highlight <- plotdat[sample==standard & same==T]
    aligned_seg$same <- aligned_seg$segment %in% same_values
    chr <- get_chr_arms()$chr
    chr <- chr[chr %in% c(1:22,'X','Y')]
    tumor_samples <- unique(aligned_seg$sample)
    tumor_samples <- tumor_samples[tumor_samples!=standard]
    cols <- c('black',rainbow(length(tumor_samples)))
    names(cols) <- c(standard, tumor_samples)
    p <- ggplot(plotdat) +
        scale_x_continuous(breaks=chr$global_midpoint, labels=chr$chr, expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        geom_vline(xintercept=chr$global_end, linewidth=0.25, color='black') +
        annotate("rect", xmin=highlight$global_seg_start, xmax=highlight$global_seg_end, ymin=min(plotdat$logR)-0.25, ymax=max(plotdat$logR)+0.25, fill='steelblue', alpha=0.08) +
        geom_segment(data=plotdat[sample!=standard], aes(x=global_seg_start, xend=global_seg_end, y=logR, yend=logR, color=sample), linewidth=0.5, alpha=0.25) +
        geom_segment(data=plotdat[sample==standard], aes(x=global_seg_start, xend=global_seg_end, y=logR, yend=logR, color=sample), linewidth=1) +
        scale_color_manual(values=cols, name='Sample') +
        facet_wrap(facets=~type, ncol=1) +
        theme_fit(base_size=10) +
        labs(x='Genomic position', y='logR')
    p    
}


# Output refphase segments as a tsv-file
write_refphase_segs <- function (segs, cn_events = NULL, file, output_format = "summary") {
    columns_to_rename <- c(group_name = "sample_id", seqnames = "chrom")
    d <- as.data.frame(segs)
    d <- S4Vectors::rename(d, columns_to_rename)
    if (output_format == "summary") {
        columns_to_save <- c("sample_id", "chrom", "start", "end",
                             "width", "cn_a", "cn_b", "was_cn_updated", "is_ai",
                             "mirrored_vs_ref", "is_reference", "any_ai", "ai_pvalue",
                             "effect_size", "diptest_pvalue", "heterozygous_SNP_number",
                             "homozygous_SNP_number")
        d_filt <- d[, columns_to_save]

    }
    else if (output_format == "full") {
        d <- merge(x = d, y = cn_events, by.x = c("sample_id",
                                                  "chrom", "start", "end", "width", "cn_a", "cn_b",
                                                  "is_ai", "mirrored_vs_ref", "is_LOH"), by.y = c("sample_id",
                                                  "seqnames", "start", "end", "width", "cn_a", "cn_b",
                                                  "is_ai", "mirrored_vs_ref", "is_LOH"))
        columns_to_save <- c("sample_id", "purity", "ploidy",
                             "chrom", "start", "end", "width", "cn_a", "cn_b",
                             "was_cn_updated", "negative_cn_called", "is_ai",
                             "mirrored_vs_ref", "is_reference", "any_ai", "ai_pvalue",
                             "effect_size", "diptest_pvalue", "heterozygous_SNP_number",
                             "homozygous_SNP_number", "amplification_logRthreshold",
                             "gain_logRthreshold", "loss_logRthreshold", "meanlogR",
                             "is_relative_amplification_mean", "is_relative_gain_mean",
                             "is_relative_loss_mean", "is_relative_amplification_ttest",
                             "is_relative_gain_ttest", "is_relative_loss_ttest",
                             "is_absolute_amplification", "is_absolute_gain",
                             "is_absolute_loss", "is_LOH", "is_homozygous_deletion")
        d_filt <- d[, columns_to_save]

    }
    else if (output_format == "copynumbers") {
        columns_to_save <- c("sample_id", "chrom", "start", "end",
                             "cn_a", "cn_b")
        d_filt <- d[, columns_to_save]
    }
    else if (output_format == "copynumbers_int") {
        columns_to_save <- c("sample_id", "chrom", "start", "end",
                             "cn_a_integer", "cn_b_integer")
        d_filt <- d[, columns_to_save]
        names(d_filt) <- gsub('_integer','',names(d_filt))
    }
    else {
        stop("output_format has to be in ['summary', 'full', 'copynumbers']")
    }
    utils::write.table(d_filt, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
}


group_cols <- c("#000000","#018B44","#EA5B2B","#F9B31C","#EA6A8B","#4B85C5","#534696","#bfbfbf")
names(group_cols) <- c('Normal','Primary','Locoregional','Peritoneum','Lung','Liver','Distant (other)','Other')


