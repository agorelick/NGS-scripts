# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here I am attempting to use pseudoSNPs + phased SNP blocks
# with mostly out of the box FACETS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
library(QDNAseq)
library(Rsamtools)
library(data.table)
library(Rcpp)
library(Biobase)
library(pctGCdata)
library(here)
library(ASCAT)

source('~/lab_repos/lpASCN/func.R')

## define input arguments/parameters
datadir <- '~/Dropbox (Partners HealthCare)/MGH CLINICAL/ALEX/germline_data/lpASCN/C161'
phased_bcf <- file.path(datadir,'original_data/C161_q=10_Q=20_P=0_slimmed_phased_merged.bcf')
pileup_data <- file.path(datadir,'original_data/C161_q=10_Q=20_P=0_slimmed.pileup.gz')
qdnaseq_data <- file.path(datadir,paste0('original_data/C161_1000kbp_withXYMT.rds'))
sample_map <- file.path(datadir,'original_data/C161_sample_map.txt') ## table with columns: pileup_order, sample_name, proper_sample_name, qdnaseq_name

d <- preprocess_data(qdnaseq_data, pileup_data, phased_bcf, sample_map)
write_tsv(d,file=file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run multipcf for multi-sample segmentation
# then extract aligned segments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))

source('~/lab_repos/lpASCN/func.R')



## below into a function that returns:
# - BAF matrix
# - logR matrix
# - bins annotated with segment number
# - wide-segments with logR and BAF
# - should it do the segmentation every time or save previous segmentation? make sure seed works also.


sweep_dipLogRs <- function(d, seg, dipLogRs, normal_sample, sex, cores) {

    sweep_sample <- function(sample, sex, seg, d, dipLogRs, cores) {
        message(sample)
        obj <- get_obj_for_sample(sample, seg, d, sex)
        sweep_dipLogR <- function(dipLogR, obj, seg, cores) {
            message(dipLogR)
            fit <- get_fit(dipLogR, obj, cores=cores)
            ascn <- get_na_nb(p=fit$pu, t=fit$pl, x=obj$sample_seg)
            ascn$dipLogR <- dipLogR
            if(!is.null(fit)) {
                ascn$len_neg <- as.numeric(fit$len_neg)
                ascn$len_homdel <-  as.numeric(fit$len_homdel)
                ascn$len_na <- as.numeric(fit$len_na)
                ascn$loglik <- as.numeric(fit$loglik)
            } else {
                ascn$len_neg <- as.numeric(NA)
                ascn$len_homdel <- as.numeric(NA)
                ascn$len_na <- as.numeric(NA)
                ascn$loglik <- as.numeric(NA)   
            }
            ascn
        }
        l <- lapply(dipLogRs, sweep_dipLogR, obj, seg, cores)
        out <- rbindlist(l)
        out$sample <- sample
        out
    }
    tumor_samples <- unique(d$sample)
    tumor_samples <- tumor_samples[tumor_samples!=normal_sample]

    l <- lapply(tumor_samples, sweep_sample, sex, seg, d, dipLogRs, cores)
    res <- rbindlist(l)
    res
}

dipLogRs <- seq(-2,2,by=0.01)
cores <- 6
this.sex <- 'female'
this.normal <- 'Normal1'

seg <- get_segments(d, normal_sample=this.normal, sex=this.sex)
res <- sweep_dipLogRs(d, seg, dipLogRs, normal_sample=this.normal, sex, cores) 


## get a matrix with total copy number for each segment in each sample/dipLogR
res$nai <- round(res$na)
res$nbi <- round(res$nb)
res[!is.na(nbi),tcn:=nai+nbi]
res[is.na(nbi),tcn:=nai]
x <- data.table::dcast(dipLogR + sample ~ segment, value.var='tcn', data=res)

plot_fit(fit, obj)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scrap below
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## why 3 clusters in BAF at chr4 segment 14?
library(Ckmeans.1d.dp)
qc <- obj$sample_dat
qc <- qc[Chromosome==4,]
qc <- qc[segment==14,]
clus <- Ckmeans.1d.dp(qc$i.BAF, k=c(1,2,3))
qc$k <- factor(clus$cluster)
p <- ggplot(qc, aes(x=bin, y=i.BAF)) +
    geom_point(aes(color=k))




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




