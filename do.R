# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here I am attempting to use pseudoSNPs + phased SNP blocks
# with mostly out of the box FACETS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
library(QDNAseq)
library(Rsamtools)
library(data.table)
library(ang)
library(Rcpp)
library(Biobase)
library(pctGCdata)


## define input arguments/parameters
datadir <- '~/Dropbox (Partners HealthCare)/MGH CLINICAL/ALEX/germline_data/lpASCN/C161'
binsize_kb <- 1000
phased_bcf <- file.path(datadir,'original_data/C161_q=10_Q=20_P=0_slimmed_phased_merged.bcf')
max_phaseable_distance <- 20000
pileup_data <- file.path(datadir,'original_data/C161_q=10_Q=20_P=0_slimmed.pileup.gz')
sample_names <- c('C161A1','C161A2','C161A5','C161A6','C161B1','C161B4','C161B6','C161B7','C161H1','C161L2','C161L3','C161L4','C161Ld1','C161Ld2','C161N1','C161O1','C161O2','C161P2','C161P3','C161P4','C161P5','C161P6','C161P7','C161P8','C161P9','C161Pa1','C161S1','C161SB1','C161SB2','C161St1','C161St2')
min_total_count_for_baf <- 5
min_allele_count_for_het <- 10


## load 1Mb bin counts from QDNAseq. We need to annotate these bins with allele-specific counts (for bins containing het-SNPs)
## this probably means the het-SNPs bins need to perfectly align with the QDNAseq bins
cnt <- readRDS(file.path(datadir,paste0('original_data/C161_',binsize_kb,'kbp_withXYMT.rds')))
cnt$sample <- gsub('_aligned.bam','',cnt$file)
cnt[,file:=NULL]
names(cnt) <- c('bin','chr','bin_start','bin_end','count','count_raw','binsize_kbp','sample')
cnt$chr <- factor(cnt$chr, levels=c(1:22,'X','Y','MT'))
cnt <- cnt[!is.na(chr),]
cnt <- cnt[order(chr,bin_start,bin_end),]

## load phased SNPs
bcf_file <- paste0(phased_bcf)
message(bcf_file)
bcf <- scanBcf(bcf_file)
phased_snps <- data.table(CHROM=bcf$CHROM, POS=bcf$POS, ID=bcf$ID, REF=bcf$REF, ALT=bcf$ALT, PHASE=bcf$GENO$GT[,1])
phased_snps <- phased_snps[CHROM %in% 1:22,]
phased_snps[,POS2:=POS+1]
phased_snps$CHROM <- factor(phased_snps$CHROM, levels=c(1:22,'X','Y','MT'))
setkey(phased_snps,'CHROM','POS','POS2')

## annotate the phased snps with 1Mb-binned counts
bins <- cnt[!duplicated(bin),c('bin','chr','bin_start','bin_end'),with=F]
setkey(bins,'chr','bin_start','bin_end')
phased_snps <- foverlaps(phased_snps, bins, type='within')
phased_snps <- phased_snps[!is.na(bin),]
phased_snps$bin <- factor(phased_snps$bin, levels=unique(phased_snps$bin))
phased_snps$bin.i <- as.integer(phased_snps$bin)

## define phasing blocks
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
get_GCbias <- function(i, phased_snps) { 
    gc <- data.table(CHROM=i, POS=phased_snps[CHROM==i,(POS)], GC=getGCpct(i, phased_snps[CHROM==i,(POS)], gbuild='hg19'))
}
gcbias <- rbindlist(lapply(1:22, get_GCbias, phased_snps))
phased_snps$GC <- gcbias$GC
phased_snps$GC <- round(phased_snps$GC, 2)

## load allele-specific read counts for the phased het-SNPs for each sample
d <- fread(pileup_data, sep=',')
d <- d[Chromosome %in% 1:22]
d$Chromosome <- factor(d$Chromosome, levels=c(1:22,'X','Y','MT'))
d <- merge(d, phased_snps[,c('CHROM','POS','PHASE','bin','block','GC'),with=F], by.x=c('Chromosome','Position'), by.y=c('CHROM','POS'), all.x=T)
d <- d[!is.na(block),]

## correct raw read count for each allele and sample for GC%
filenames <- grep('File',names(d),value=T)
for(f in filenames) {
    message(f)
    x <- d[,c('GC',filename),with=F]
    names(x)[2] <- 'count'
    qs <- quantile(x$count,c(0.005,0.995),na.rm=T)
    x[(count < qs[1] | count > qs[2]), count:=NA]
    grand_mean <- mean(x$count, na.rm=T)
    x[,diff_from_mean:=count - grand_mean]
    fit <- lm(diff_from_mean ~ GC, data=x)
    predicted_diff_from_mean <- predict(fit, newdata=d, type='response')
    d[[f]] <- d[[f]]-predicted_diff_from_mean
}

## replace File-names with sample-names
file_fields <- grep('File',names(d),value=T)
files <- unique(substr(file_fields, 1, nchar(file_fields)-1))
filesA <- data.table(newname=paste0(sample_names,'A'), file=paste0(files,'A'))
filesR <- data.table(newname=paste0(sample_names,'R'), file=paste0(files,'R'))
files <- rbind(filesA, filesR)
d_fields <- data.table(pos=1:ncol(d), name=names(d))
d_fields <- merge(d_fields, files, by.x='name', by.y='file', all.x=T)
d_fields[is.na(newname), newname:=name]
d_fields <- d_fields[order(pos),]
names(d) <- d_fields$newname

## extract the corrected read counts per block for each allele
count10_A <- d[PHASE=='1|0',c('block',paste0(sample_names,'A')),with=F]; names(count10_A) <- c('block',sample_names)
count01_R <- d[PHASE=='0|1',c('block',paste0(sample_names,'R')),with=F]; names(count01_R) <- c('block',sample_names)
count1 <- rbind(count10_A, count01_R)
count01_A <- d[PHASE=='0|1',c('block',paste0(sample_names,'A')),with=F]; names(count01_A) <- c('block',sample_names)
count10_R <- d[PHASE=='1|0',c('block',paste0(sample_names,'R')),with=F]; names(count10_R) <- c('block',sample_names)
count2 <- rbind(count01_A, count10_R)
collapse_block <- function(count) {
    out <- colSums(count)     
    as.list(out)
}
block_counts1 <- count1[,collapse_block(.SD),by=c('block')]
block_counts2 <- count2[,collapse_block(.SD),by=c('block')]


## generate the BAF matrix
block_counts <- merge(block_counts1, block_counts2, by='block', all=F) ## merge the corrected counts for each allele per block. Only retain blocks with counts for both alleles
baf <- matrix(nrow=nrow(block_counts), ncol=length(sample_names))
rownames(baf) <- block_counts$block
colnames(baf) <- sample_names
for(s in sample_names) {
    message(s)
    baf[,s] <- block_counts[[paste0(s,'.x')]] / (block_counts[[paste0(s,'.x')]] + block_counts[[paste0(s,'.y')]])
}
bad_blocks <- which(is.na(rowSums(baf)))
baf <- baf[-bad_blocks, ]
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
d <- d[,c('chr','position','Ref','Alt','sample','count1','count2','count_raw'),with=F]
d <- d[order(chr,position,sample),]


## update the sample names and save this as prepped data
sample_info <- fread('~/lab_repos/crc_lung_met/processed_data/sample_info.txt')
sample_info <- sample_info[Patient_ID=='C161',c('Sample_ID','Real_Sample_ID'),with=F]
d <- merge(sample_info, d, by.x='Sample_ID', by.y='sample', all.y=T)
write_tsv(d,file=file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_prepped.tsv'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run multipcf for multi-sample segmentation
# then extract aligned segments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ASCAT)
source('~/lab_repos/lpASCN/run_multipcf.R')

## load the prepped data and run ASCAT's multipcf segmentation algorithm
d <- fread(file=file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_prepped.tsv'))
d[is.na(count_raw),count_raw:=0]
d[,dp:=count1+count2]
d[Ref=='N', BAF:=count1 / dp]
d[,bin:=paste0(chr,':',position)]
d$bin <- as.integer(factor(d$bin, levels=unique(d$bin)))
get_logR <- function(d) {
    mid <- median(d$dp,na.rm=T)
    d$logR <- log2(d$dp / mid)
    d
}
d <- d[,get_logR(.SD),by=Real_Sample_ID]
setnames(d,c('chr','position','Real_Sample_ID'),c('Chromosome','Position','sample'))
samples <- unique(d$sample)
samples <- samples[samples!='Normal1']
run_multipcf(d, sex='female', build='hg19', penalty=70, tmpdir=file.path(datadir,'processed_data/facets/multipcf'), cleanup=F, seed=42, output_dir=file.path(datadir,'processed_data/facets/multipcf'))
## this is capping the logR in chrM, need to fix this.


## load the multipcf data and get logR and BAF matrices
load_multipcf_data <- function(sample) {
    message(sample)
    multipcf_file <- file.path(datadir,paste0('processed_data/facets/multipcf/',sample,'_multipcf.txt'))
    pcf <- fread(multipcf_file)
    pcf
}
l <- rbindlist(lapply(samples, load_multipcf_data))
l[,bin:=bin+1]
l <- l[Chromosome %in% c(1:22,'X','Y','MT'),]
l$Chromosome <- factor(l$Chromosome, levels=c(1:22,'X','Y','MT'))
l <- l[order(Chromosome, Position, sample),]
logR <- data.table::dcast(bin + Chromosome + Position ~ sample, value.var='logR', data=l)
bins <- logR[,c('bin','Chromosome','Position'),with=F]
logR[,c('Chromosome','Position'):=NULL]
logR <- ang::d2m(logR)
BAF <- ang::d2m(data.table::dcast(bin ~ sample, value.var='BAF', data=l))


## get aligned segments across all samples
n <- nrow(BAF)
this.logR <- logR[1,]
this.BAF <- BAF[1,]
segment <- rep(as.integer(NA), n)
segment[1] <- 1
segment.count <- 1
for(i in 2:n) {
    next.logR <- logR[i,]
    next.BAF <- BAF[i,]
    valid.logR <- !is.na(next.logR) & !is.na(this.logR)
    valid.BAF <- !is.na(next.BAF) & !is.na(this.BAF)
    if( any(next.logR[valid.logR]!=this.logR[valid.logR]) | any(next.BAF[valid.BAF]!=this.BAF[valid.BAF]) ) {
        segment.count <- segment.count + 1
    }
    this.logR <- next.logR; this.BAF <- next.BAF
    segment[i] <- segment.count + 1 
}
bins$segment <- segment


## get wide segment-level data for logR and BAF
ll <- merge(l, bins[,c('bin','segment'),with=F], all.x=T)
collapse_segment <- function(ll) {
    seg_start <- min(ll$Position)
    seg_end <- max(ll$Position)
    BAF <- unique(ll$BAF[!is.na(ll$BAF)])
    logR <- unique(ll$logR[!is.na(ll$logR)])
    if(length(BAF) > 1 | length(BAF) > 1) stop('Multiple values per segment?')
    list(seg_start=seg_start, seg_end=seg_end, BAF=BAF, logR=logR)
}
seg <- ll[,collapse_segment(.SD), by=c('sample','Chromosome','segment')]


source('~/lab_repos/lpASCN/get_chr_arms.R')
chr_arms <- get_chr_arms()
chr <- chr_arms$chr
seg <- merge(seg, chr[,c('chr','global_start'),with=F], by.x='Chromosome', by.y='chr', all.x=T)
seg[,global_seg_start:=seg_start/1e6+global_start]
seg[,global_seg_end:=seg_end/1e6+global_start]
seg[,seg_length:=seg_end - seg_start + 1]



get_na_nb <- function(p, t, x) {
    ## I think this works?
    r <- x$logR
    b <- x$BAF
    na = -1*((b*p*t*2^r) - (b*p*2^(r+1)) + (b*2^(r+1)) - (p*t*2^r) + (p*2^(r+1)) - p - 2^(r+1) + 1) / p
    nb = ((b*p*t*2^r) - (b*p*2^(r+1)) + (b*2^(r+1)) + p - 1) / p

    ## go over previous values where BAF was NA
    NAs <- is.na(b)
    c <- 2^x$logR[NAs]
    na[NAs] <- (2*(c-1) / p) + t*c - 2*(c-1)
    nb[NAs] <- as.numeric(NA)

    x$na <- na
    x$nb <- nb
    x
}


get_logR <- function(pu, pl, na, nb) log2((2*(1-pu) + pu*(na + nb)) / (2*(1-pu) + pu*pl))
get_BAF <- function(pu, pl, na, nb) (1 - pu + pu*nb) / (2 - 2*pu + pu*(na+nb))


get_loglik <- function(pu, pl, sample_seg, sample_dat) {
    sample_seg2 <- get_na_nb(p=pu, t=pl, sample_seg)
    sample_seg2 <- sample_seg2[!is.na(na) & !is.na(nb) & !is.na(BAF) & !is.na(logR),]
    sample_seg2$na <- round(sample_seg2$na)
    sample_seg2$nb <- round(sample_seg2$nb)
    sample_seg2[na < 0, na:=0]
    sample_seg2[nb < 0, nb:=0]
    sample_seg2$logR_expected <- get_logR(pu, pl, sample_seg2$na, sample_seg2$nb)
    sample_seg2$BAF_expected <- get_BAF(pu, pl, sample_seg2$na, sample_seg2$nb)
    sample_dat2 <- merge(sample_dat, sample_seg2[,c('segment','logR_expected','BAF_expected'),with=F], by='segment', all.x=T)
    sample_dat2 <- sample_dat2[!is.na(logR_expected) & !is.na(BAF_expected),] 

    ## get the expected logR and BAF values if na,nb are the nearest integers
    diff_from_logR_expected <- sample_dat2$i.logR - sample_dat2$logR_expected
    diff_from_logR_expected <- diff_from_logR_expected[!is.na(diff_from_logR_expected)]
    diff_from_BAF_expected <- sample_dat2$i.BAF - sample_dat2$BAF_expected
    diff_from_BAF_expected <- diff_from_BAF_expected[!is.na(diff_from_BAF_expected)]

    ## calculate the summed loglikelihood of the bins' observed logRs and BAFs, given the nearest integer copy numbers
    loglik_logR <- log(c(logR_CDF(diff_from_logR_expected[diff_from_logR_expected >= 0]), logR_CDF(diff_from_logR_expected[diff_from_logR_expected < 0])))
    min_valid_logR_loglik <- min(loglik_logR[!is.infinite(loglik_logR)])
    loglik_logR[is.infinite(loglik_logR)] <- min_valid_logR_loglik
    loglik_BAF <- log(c(BAF_CDF(diff_from_BAF_expected[diff_from_BAF_expected >= 0]), BAF_CDF(diff_from_BAF_expected[diff_from_BAF_expected < 0])))
    min_valid_BAF_loglik <- min(loglik_BAF[!is.infinite(loglik_BAF)])
    loglik_BAF[is.infinite(loglik_BAF)] <- min_valid_BAF_loglik

    loglik <- sum(loglik_logR) + sum(loglik_BAF)
    loglik
}


## plot the segments with the raw logR/BAF values to see what I'm working with.
get_fit <- function(dipLogR, sample_dat, sample_seg, max_homdel=1e8, max_neg=0, max_ploidy=8) {
    #max_homdel=1e8; max_neg=0; max_ploidy=8; dipLogR <- -0.432

    require(parallel)
    purity <- seq(0.05,1,by=0.001)
    ploidy <- (2*(purity + 2^(-dipLogR) - 1)) / purity
    myfits <- data.table(pu=purity, pl=round(ploidy,3))

    ## for all fits, first filter out impossible fits based on negative copies and too much homozygous deletion
    qc_fit <- function(i, myfits, sample_seg) {
        this.pu <- myfits$pu[i]
        this.pl <- myfits$pl[i]
        fit <- get_na_nb(this.pu, this.pl, sample_seg)
        len_neg <- sum(sample_seg$seg_length[round(fit$na) <= -1 | round(fit$nb) <= -1],na.rm=T)
        len_homdel <- sum(sample_seg$seg_length[round(fit$na) == 0 & round(fit$nb) == 0],na.rm=T)
        len_na <- sum(sample_seg$seg_length[is.na(fit$na) | is.na(fit$nb)])
        list(pu=this.pu, pl=this.pl, len_neg=len_neg, len_homdel=len_homdel, len_na=len_na)
    }

    n_fits <- nrow(myfits)
    profile_list <- mclapply(1:n_fits, qc_fit, myfits, sample_seg, mc.cores=4) 
    myfits <- rbindlist(profile_list)
    myfits <- myfits[len_neg <= max_neg & len_homdel <= max_homdel & pl <= max_ploidy,]

    get_possible_fits_for_sample <- function(myfits, sample_seg, sample_dat, cores=4) {
        try_fits <- function(i, myfits, sample_seg, sample_dat) {
            pu <- myfits$pu[i]
            pl <- myfits$pl[i]
            loglik <- get_loglik(pu, pl, sample_seg, sample_dat)
            list(pl=pl, pu=pu, loglik=loglik)
        }
        n_fits <- nrow(myfits)
        l <- mclapply(1:n_fits, try_fits, myfits, sample_seg, sample_dat, mc.cores=cores)
        res <- rbindlist(l)
        myfits$loglik <- res$loglik
        myfits
    }

    if(nrow(myfits) > 0) {
        myfits <- get_possible_fits_for_sample(myfits, sample_seg, sample_dat)
        myfits <- myfits[order(loglik,decreasing=T),]
        myfits$dipLogR <- dipLogR
        myfits[1,]
    } else {
        NULL
    }
}


this.sample <- 'PT8-A'

sample_seg <- seg[sample==this.sample]
sample_seg[,seg_end:=seg_end+(500000-1)]
sample_seg[,seg_start:=seg_start-(500000-1)]
sample_dat <- d[sample==this.sample]
sample_dat[,Position2:=Position+1]
setkey(sample_dat,'Chromosome','Position','Position2')
toadd <- sample_seg[,c('Chromosome','seg_start','seg_end','BAF','logR','segment','global_start'),with=F]
setkey(toadd,'Chromosome','seg_start','seg_end')
sample_dat <- foverlaps(sample_dat, toadd, type='within')
sample_dat$Chromosome <- factor(sample_dat$Chromosome, levels=c(1:22,'X','Y','MT'))
sample_dat <- sample_dat[order(Chromosome, Position),]
sample_dat[,global_pos:=Position/1e6+global_start]
sample_dat[,c('seg_start','seg_end'):=NULL]

## in each segment, what does the distribution around the logR and BAF estimates look like?
qc <- sample_dat[,c('segment','i.BAF','i.logR','BAF','logR'),with=F]
qc[,diff_from_BAF:=i.BAF - BAF]
qc[,diff_from_logR:=i.logR - logR]

#library(EnvStats)
#demp(qc$diff_from_logR)
logR_CDF <- ecdf(qc$diff_from_logR)
BAF_CDF <- ecdf(qc$diff_from_BAF)



#fit <- get_fit(-0.73, sample_dat, sample_seg, max_homdel=1e8, max_neg=0) 
#fit <- get_fit(-1, sample_dat, sample_seg, max_homdel=1e8, max_neg=0) 
#fit <- get_fit(0.3, sample_dat, sample_seg, max_homdel=1e8, max_neg=0) 


## currently this is returning -Inf loglik's because there can be more extreme diff-from-expected-BAF than empirically observed.
## how do I prevent them from having prob=0?
fit <- get_fit(-0.432, sample_dat, sample_seg, max_homdel=1e8, max_neg=0) 


plot_fit <- function(fit, sample_seg, sample_dat) {
    require(cowplot)
    require(ggplot2)

    plot_sample_seg <- sample_seg[Chromosome %in% c(1:22,'X')]
    plot_sample_dat <- sample_dat[Chromosome %in% c(1:22,'X')]
    plot_chr <- chr[chr %in% c(1:22,'X')]
    plot_chr$chr <- factor(plot_chr$chr, levels=c(1:22,'X'))
    plot_sample_dat$chr <- factor(plot_sample_dat$chr, levels=c(1:22,'X'))
    plot_sample_seg$chr <- factor(plot_sample_seg$chr, levels=c(1:22,'X'))
    centerline <- median(plot_sample_seg$logR,na.rm=T)
    dipLogR <- fit$dipLogR

    p1 <- ggplot(plot_sample_dat, aes(x=global_pos, y=i.logR)) + 
        geom_point(color='#bfbfbf',size=1.25,pch=16) +
        geom_hline(yintercept=centerline, color='green', size=0.25) +
        geom_hline(yintercept=dipLogR, color='purple', size=0.25) +
        geom_vline(xintercept=c(0,plot_chr$global_end)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr, expand=c(0,0)) +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start,xend=global_seg_end,y=logR,yend=logR),color='blue',linewidth=0.75) +
        theme_ang(base_size=12)
    p2 <- ggplot(plot_sample_dat, aes(x=global_pos, y=i.BAF)) + 
        geom_point(color='#bfbfbf',size=1.25,pch=16) +
        geom_vline(xintercept=c(0,plot_chr$global_end)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr,expand=c(0,0)) +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start,xend=global_seg_end,y=BAF,yend=BAF),color='blue',linewidth=0.75) +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start,xend=global_seg_end,y=1-BAF,yend=1-BAF),color='blue',linewidth=0.75) +
        theme_ang(base_size=12)
    ascn <- get_na_nb(p=fit$pu, t=fit$pl, x=plot_sample_seg)
    ascn[,tcn:=na+nb]
    ascn[na >= nb,mcn:=nb]
    ascn[na < nb,mcn:=na]
    ascn[tcn >= 10, tcn:=log10(tcn)+9]
    ascn[mcn >= 10, tcn:=log10(mcn)+9]
    ascn[tcn < -1, tcn:=-1]
    ascn[mcn < -1, mcn:=-1]

    p3 <- ggplot(ascn) + 
        scale_y_continuous(breaks=c(seq(0,10,by=2),log10(50)+9,log10(100)+9), labels=c(seq(0,10,by=2),50,100)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr, expand=c(0,0)) +
        geom_hline(yintercept=seq(0,10,by=1), color='#bfbfbf', size=0.25) +
        geom_vline(xintercept=c(0,plot_chr$global_end),size=0.25) +
        geom_segment(aes(x=global_seg_start,xend=global_seg_end,y=tcn,yend=tcn),color='black',linewidth=0.75) +
        geom_segment(aes(x=global_seg_start,xend=global_seg_end,y=mcn,yend=mcn),color='red',linewidth=0.75) +
        theme_ang(base_size=12)

    p <- plot_grid(p1, p2, p3, align='v', ncol=1)
    p
}

fit <- get_fit(-0.43, sample_dat, sample_seg, max_homdel=1e8, max_neg=0) 
plot_fit(fit, sample_seg, sample_dat)
ascn <- get_na_nb(p=fit$pu, t=fit$pl, x=sample_seg)

#sample_fits <- fits[sample=='PT8-A']
















get_possible_fits_for_sample <- function(this.sample, seg, cores=4) {
    require(parallel)
    message(this.sample)
    purity <- round(seq(0.05,0.95,by=0.01),2)
    ploidy <- round(seq(1.0, 6, by=0.01),2)
    pupl <- as.data.table(expand.grid(pu=purity, pl=ploidy))
    n_fits <- nrow(pupl)

    try_fits <- function(i, pupl, x) {
        pu <- pupl$pu[i]
        pl <- pupl$pl[i]
        MSE <- test_pupl(pu, pl, x)
        list(pl=pl, pu=pu, MSE=MSE)
    }

    x <- seg[sample==this.sample]
    l <- mclapply(1:n_fits, try_fits, pupl, x, mc.cores=cores)
    res <- rbindlist(l)
    res$sample <- this.sample
    res
}
fit_list <- lapply(samples, get_possible_fits_for_sample, seg)
fits <- rbindlist(fit_list)
fits$dipLogR <- log2(2/(2*(1-fits$pu) + fits$pu*fits$pl))
fits[sample=='PT8-A']

## remove bad fits based on negative copy number segments
get_profile_for_fit <- function(i, fits, seg) {
    this.sample <- fits$sample[i]
    this.pu <- fits$pu[i]
    this.pl <- fits$pl[i]
    this.dipLogR <- fits$dipLogR[i]
    sample.seg <- seg[sample==this.sample]
    fit <- get_na_nb(this.pu, this.pl, sample.seg)
    total_length <-  sum(sample.seg$seg_length)
    len_negative <- sum(sample.seg$seg_length[round(fit$na) <= -1 | round(fit$nb) <= -1],na.rm=T)
    len_homdel <- sum(sample.seg$seg_length[round(fit$na) == 0 & round(fit$nb) == 0],na.rm=T)
    len_na <- sum(sample.seg$seg_length[is.na(fit$na) | is.na(fit$nb)])
    list(sample=this.sample, purity=this.pu, ploidy=this.pl, dipLogR=this.dipLogR, len_negative=len_negative, len_homdel=len_homdel, len_na=len_na, total_length=total_length)    
}
n_fits <- nrow(fits)
#n_fits <- 1e3
profile_list <- mclapply(1:n_fits, get_profile_for_fit, fits, seg, mc.cores=4) ## ~30 sec
profiles <- rbindlist(profile_list)






get_fit_for_dipLogR <- function(dipLogR, sample_fits) {
    sample_fits$diff_from_dipLogR <- abs(sample_fits$dipLogR - dipLogR)
    sample_fits <- sample_fits[order(diff_from_dipLogR, decreasing=F),]
    sample_fits <- sample_fits[diff_from_dipLogR==min(diff_from_dipLogR),]
    sample_fits <- sample_fits[MSE==min(MSE),]
    sample_fits
}
myfit <- get_fit_for_dipLogR(dipLogR, sample_fits)












## for each ploidy value, get the most-likely purity estimates based on local minima of MSE
get_local_minima <- function(x) {
    minima <- which(diff(sign(diff(x$MSE)))==2)+1
    x[minima]
}
minima_fits <- fits[,get_local_minima(.SD),by=c('sample','pl')]









## for each fit, check the number of segments that h
#samples %in% profiles$sample ## we still have all samples

## for each sample, get the lower quantile of the MSE, use this to reduce the number of minima fits
#sample_thresholds <- function(fits) {
#    q=quantile(fits$MSE,0.2)
#    list(q=q)
#}
#qs <- fits[,sample_thresholds(.SD),by=sample]
#minima_fits <- merge(minima_fits, qs, by='sample', all.x=T)
#minima_fits <- minima_fits[MSE < q,]

minima_fits$minimum <- T
pd <- merge(fits, minima_fits[,c('sample','pl','pu','minimum'),with=F], by=c('sample','pl','pu'), all.x=T)
pd[is.na(minimum),minimum:=F]
p <- ggplot(pd, aes(x=pl, y=pu)) +
    geom_tile(aes(fill=MSE)) +
    geom_point(data=pd[minimum==T,], pch=16, size=0.1, color='black') +
    facet_wrap(facets=~sample) +
    theme_ang(base_size=12) +
    scale_fill_gradient2(low='red',mid='white',high='blue',midpoint=median(minima_fits$MSE))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 









x <- seg[sample=='PT8-A']
#x[,logR:=logR+0.60]
#x[,logR:=logR+0.54]
#x[,logR:=logR-0.36]
ascn <- get_na_nb(p=0.44, t=3.3, x)
x$na <- ascn$na
x$nb <- ascn$nb
#x$na <- round(ascn$na)
#x$nb <- round(ascn$nb)
x[Chromosome==3,]
x[,n:=round(na)+round(nb)]
sum(x$n*x$seg_length) / sum(x$seg_length)
#profiles[sample=='PT8-A' & ploidy==3.3 & purity==0.73]





get_profile_for_fit <- function(i, minima_fits, seg) {
    this.sample <- minima_fits$sample[i]
    this.pu <- minima_fits$pu[i]
    this.pl <- minima_fits$pl[i]
    sample.seg <- seg[sample==this.sample]
    fit <- get_na_nb(this.pu, this.pl, sample.seg)
    fit$purity <- this.pu
    fit$ploidy <- this.pl
    fit$nai <- round(fit$na)
    fit$nbi <- round(fit$nb)
    out <- as.list(c(fit$nai, fit$nbi))
    names(out) <- c(paste0(fit$segment,'.a'), paste0(fit$segment,'.b'))
    out <- append(list(sample=this.sample, purity=this.pu, ploidy=this.pl), out)
    out
}
profile_list <- mclapply(1:nrow(minima_fits), get_profile_for_fit, minima_fits, seg, mc.cores=4) ## ~30 sec
profiles <- rbindlist(profile_list)

## try UMAP to see how samples cluster based on integer copy number (among the positive-minima-fits)
set.seed(42)
library(umap)
sample_profiles <- profiles
umap.data = as.matrix(sample_profiles[,(4:ncol(sample_profiles)),with=F])
bad_cols <- which(is.na(colSums(umap.data)))
umap.data <- umap.data[,-bad_cols]
my.umap = umap(umap.data, n_components = 2)
layout <- my.umap[["layout"]] 
layout <- as.data.table(layout) 
final <- cbind(layout, sample_profiles[,c(1:3),with=F])
best <- final[(V1 > -6 & V1 < 1) & (V2 > -1 & V2 < 6)] ## this gives many samples with ploidy ~ 1
p <- ggplot(final, aes(x=V2, y=V1)) +
    geom_point(aes(fill=sample),pch=21,color='black',stroke=0.25) 


    facet_wrap(facets=~sample)

## can I try to identify segments to target with the linear program??
dm <- data.table::melt(profiles, id.vars=c('sample','purity','ploidy'))
dm$variable <- as.character(dm$variable)
dm[,segment:=as.integer(strtrim(variable,(nchar(variable)-2)))]
dm[,allele:=substr(variable,(nchar(variable)),(nchar(variable)))]
dm <- data.table::dcast(sample + purity + ploidy + segment ~ allele, value.var='value', data=dm)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For comparison, run FACETS for samples individually
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

binsize_kb <- 1000
max_phaseable_distance <- 20000
sample_names <- c('C161A1','C161A2','C161A5','C161A6','C161B1','C161B4','C161B6','C161B7','C161H1','C161L2','C161L3','C161L4','C161Ld1','C161Ld2','C161N1','C161O1','C161O2','C161P2','C161P3','C161P4','C161P5','C161P6','C161P7','C161P8','C161P9','C161Pa1','C161S1','C161SB1','C161SB2','C161St1','C161St2')
min_total_count_for_baf <- 5
min_allele_count_for_het <- 10
source('~/lab_repos/facets/multiplot.R')

dat <- fread(file='processed_data/C161_q=10_Q=20_P=0_1000kbs_prepped.tsv')
setnames(dat,c('chr','position'),c('Chromosome','Position'))
dat[,sample:=Real_Sample_ID]
dat <- dat[Chromosome %in% c(1:22,'X')]

## format into snp-bins
dat$Chromosome <- factor(dat$Chromosome, levels=c(1:22,'X'))

## pull out the positions that are heterozygous
het.pos <- as.data.frame(dat[sample=='Normal1',c('Chromosome','Position','count1','count2','Ref'),with=F])
het.pos$het <- as.integer(het.pos$Ref=='N' & het.pos$count1 > min_allele_count_for_het & het.pos$count2 > min_allele_count_for_het)
het.pos <- het.pos[het.pos$Chromosome %in% 1:22,c('Chromosome','Position','het')]
het.pos$Chromosome <- as.integer(het.pos$Chromosome)

Ndat <- dat[sample=='Normal1',c('Chromosome','Position','Ref','Alt','count1','count2'),with=F]
setnames(Ndat,c('count1','count2'),c('File1R','File1A'))
Ndat$File1E <- 0
Ndat$File1D <- 0

get_facets_input_for_sample <- function(this.sample, d) {
    message(this.sample)
    Tdat <- dat[sample==this.sample,c('Chromosome','Position','Ref','Alt','count1','count2'),with=F]
    setnames(Tdat,c('count1','count2'),c('File2R','File2A'))
    Tdat$File2E <- 0
    Tdat$File2D <- 0

    ## merge the two samples' data, only include the bins which had values in both the current tumor and the normal
    out <- merge(Ndat, Tdat, by=c('Chromosome','Position','Ref','Alt'), all=F)
    out[is.na(File1R),File1R:=0]
    out[is.na(File1A),File1A:=0]
    out[is.na(File2R),File2R:=0]
    out[is.na(File2A),File2A:=0]
    out <- out[order(Chromosome, Position),]
    patient_dir <- 'processed_data/facets'
    if(!dir.exists(patient_dir)) dir.create(patient_dir, recursive=T)
    write.table(out, file = paste0('processed_data/facets/inputdata_',this.sample,'_facets_input.csv'), sep = ',', quote = F, row.names = F)
}
samples <- sort(unique(dat[sample!='Normal1',(sample)]))
trash <- lapply(samples, get_facets_input_for_sample, d)


get_fit <- function(this.sample, cval, dipLogR=NA, save=F, seed=1234, return.obj=F, facets.gc.correction=F, cf.em.ymax=10, het.pos=NULL, max_abs_logOR=5) {

    #this.sample <- 'PT8-A'; cval <- 50; facets.gc.correction <- F
    inputdata <- paste0('processed_data/facets/inputdata_',this.sample,'_facets_input.csv')
    rcmat = readSnpMatrix(inputdata)
    rcmat <- rcmat[rcmat$Chromosome %in% 1:22,]
    rcmat$Chromosome <- as.integer(rcmat$Chromosome)
    #rcmat$Chromosome <- factor(rcmat$Chromosome, levels=c(1:22,'X'))
    rcmat <- rcmat[order(rcmat$Chromosome, rcmat$Position),]
    rcmat <- rcmat[!is.na(rcmat$NOR.DP) & !is.na(rcmat$TUM.DP),]   

    ## intercept the input data for facets to remove problematic values
    qc <- as.data.table(rcmat)
    qc[,logOR:=log2( (TUM.RD / (TUM.DP-TUM.RD)) / (NOR.RD / (NOR.DP-NOR.RD)) )]
    qc <- qc[abs(logOR) <= max_abs_logOR | is.nan(logOR),]
    rcmat <- as.data.frame(qc[,c(1:6),with=F])
    xx = preProcSample(rcmat, ndepth=2*min_allele_count_for_het, ndepthmax=Inf, het.pos=het.pos, correct.gc=F, hetscale=T)

    if(is.na(dipLogR)) {
        oo=procSample(xx,cval=cval)
        list(dipLogR=oo$dipLogR, alBalLogR=oo$alBalLogR, flags=oo$flags)

    } else if(!is.na(dipLogR) & save==F) {
        oo=procSample(xx,cval=cval, dipLogR=dipLogR)
        fit=emcncf(oo)
        fit$purity <- round(fit$purity,3)
        fit$ploidy <- round(fit$ploidy,3)
        fit$dipLogR <- round(fit$dipLogR,3)
        fit$loglik <- round(fit$loglik,3)
        message('purity: ',fit$purity)
        message('ploidy: ',fit$ploidy)
        message('dipLogR: ',oo$dipLogR)
        sname <- paste0(this.sample,' purity=',round(fit$purity,3),'; ploidy=',round(fit$ploidy,3),'; dipLogR=',round(oo$dipLogR,3))
        if(return.obj==T) {
            out <- merge(oo$jointseg, fit$cncf[,c('seg','tcn.em','lcn.em','cf.em')], by='seg', all.x=T)
            out$sample <- this.sample
            out <- out[,c('sample','seg','chrom','maploc','tcn.em','lcn.em','cf.em')]
            return(out) 
        } else {
            multiplot(x=oo,emfit=fit,sname=sname, cf.em.ymax=cf.em.ymax, tcn.field='tcn.em', lcn.field='lcn.em') 
        }
    } else if(!is.na(dipLogR) & save==T) {
        oo=procSample(xx,cval=cval, dipLogR=dipLogR)
        fit=emcncf(oo)
        sname <- paste0(this.sample,' purity=',round(fit$purity,3),'; ploidy=',round(fit$ploidy,3),'; dipLogR=',round(oo$dipLogR,3))

        ## get file names
        fit_file <- paste0('processed_data/facets/fit_',this.sample,'_purity=',round(fit$purity,3),'_ploidy=',round(fit$ploidy,3),'_dipLogR=',round(fit$dipLogR,3),'.pdf')
        obj_file <- paste0('processed_data/facets/rds_',this.sample,'_purity=',round(fit$purity,3),'_ploidy=',round(fit$ploidy,3),'_dipLogR=',round(fit$dipLogR,3),'.rds')

        pdf(fit_file, height=10.5, width=8)
        multiplot(x=oo,emfit=fit,sname=sname, cf.em.ymax=cf.em.ymax, tcn.field='tcn.em', lcn.field='lcn.em')
        dev.off()

        # Generate output
        obj = list(
                   tumor_sample = this.sample,
                   normal_sample = 'Normal1',
                   snps = oo$jointseg,
                   segs = fit$cncf,
                   out = oo$out,
                   purity = as.numeric(fit$purity),
                   ploidy = as.numeric(fit$ploidy),
                   dipLogR = oo$dipLogR,
                   alBalLogR = oo$alBalLogR,
                   flags = oo$flags,
                   em_flags = fit$emflags,
                   loglik = fit$loglik,
                   seed=seed,
                   raw_input = rcmat
        )
        saveRDS(obj, file=obj_file)
        if(return.obj==T) return(obj)
    } else {
        message('bad entry.') 
    }
}

library(facets)
this.sample <- 'PT8-A'
obj <- get_fit(this.sample, cval=100, dipLogR=NA, save=F, seed=1234, return.obj=F, facets.gc.correction=F, het.pos=het.pos)
pdf('~/Dropbox (Partners HealthCare)/Naxerova lab/Alex/Presentations/20231002_lab_meeting_lpASCN/PT8-A_facets_fit1.pdf')
get_fit(this.sample, cval=25, dipLogR=-0.3946854, save=F, seed=1234, return.obj=F, facets.gc.correction=F, het.pos=het.pos)
dev.off()

#get_fit(this.sample, cval=100, dipLogR=-0.54, save=F, seed=1234, return.obj=F, facets.gc.correction=F, het.pos=het.pos)
#dev.off()

x <- fits[sample==this.sample,]
x$pl2 <- round(x$pl,1)
p <- ggplot(x[pu==2.8,], aes(x=pu, y=MSE)) + 
    geom_line(color='steelblue') +
    geom_point(color='steelblue',pch=16,size=1.5) +
    #geom_vline(xintercept=0.718,color='red',size=0.5,linetype='dashed') +
    theme_ang(base_size=12) +
    labs(x='Purity', y='Mean squared-error', title='PT8-A purity, given ploidy=2.8')
ggsave('~/Dropbox (Partners HealthCare)/Naxerova lab/Alex/Presentations/20231002_lab_meeting_lpASCN/PT8-A_mse.pdf',width=8,height=6)


fits[sample==this.sample & pu==0.74]
fits[sample==this.sample & pu >= 0.6,(pl)]







autofit <- function(sample) {
    message(sample)
    obj <- get_fit(sample, cval=100, dipLogR=NA, save=F, seed=1234, return.obj=F, facets.gc.correction=F, het.pos=het.pos)
    obj <- get_fit(sample, cval=100, dipLogR=obj$dipLogR, save=T, seed=1234, return.obj=T, facets.gc.correction=F, het.pos=het.pos)
    obj
}






l <- lapply(samples, autofit)
saveRDS(l,file='processed_data/facets/all_fits.rds')

#pdf(file=paste0('testfit_C161_PT8-A_binsize=',binsize_kb,'kb_dipLogR=-0.52_cval=1000.pdf'))
#get_fit('PT8-A', cval=1000, dipLogR=-0.52, save=F, seed=1234, return.obj=F, facets.gc.correction=F, cf.em.ymax=NA, het.pos=het.pos)
#dev.off()
#pdf(file=paste0('testfit_C161_PT8-A_binsize=',binsize_kb,'kb_dipLogR=-0.52_cval=500.pdf'))
#get_fit('PT8-A', cval=500, dipLogR=-0.52, save=F, seed=1234, return.obj=F, facets.gc.correction=F, cf.em.ymax=NA, het.pos=het.pos)
#dev.off()
#pdf(file=paste0('testfit_C161_PT8-A_binsize=',binsize_kb,'kb_dipLogR=-0.52_cval=100.pdf'))
#get_fit('PT8-A', cval=100, dipLogR=-0.52, save=F, seed=1234, return.obj=F, facets.gc.correction=F, cf.em.ymax=NA, het.pos=het.pos)
#dev.off()




