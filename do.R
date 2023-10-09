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
seg <- get_segments(d, normal_sample='Normal1', sex='female')
write_tsv(seg,file=file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_segments.tsv'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run multipcf for multi-sample segmentation
# then extract aligned segments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


source('~/lab_repos/lpASCN/func.R')
d <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))

seg <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_segments.tsv'))
dat <- d2m(data.table::dcast(segment ~ sample, value.var='logR', data=seg))
write_tsv(dat,file='/Users/alexgorelick/lab_repos/lpASCN/C161_data_for_gurobi.tsv')

## choose our 'standard' sample
standard <- 'PT8-A'
valid_segs <- seg[seg_length >= 1e7 & sample==standard] ## subset 
valid_segs <- valid_segs[order(seg_length,decreasing=T),]
valid_segs <- valid_segs[!duplicated(Chromosome),]
seg <- seg[segment %in% valid_segs$segment,]

## first get the threshold error for a segment to be considered to have the 'same' logR in each sample
tmp.d <- d[,c('sample','bin','logR','Chromosome','Position'),with=F]
tmp.d[,Position2:=Position+1]
tmp.d[,id:=paste0(sample,'.',Chromosome)]
seg.d <- seg[,c('sample','segment','Chromosome','seg_start','seg_end','logR'),with=F]
seg.d[,id:=paste0(sample,'.',Chromosome)]
setkey(tmp.d,'id','Position','Position2')
setkey(seg.d,'id','seg_start','seg_end')
tmp.d <- foverlaps(tmp.d, seg.d, type='within')
tmp.d <- tmp.d[!is.na(segment),]
tmp.d <- tmp.d[Chromosome %in% c(1:22,'X')]
tmp.d[,diff_from_seg:=i.logR - logR]
quantile(abs(tmp.d$diff_from_seg),0.20)

## export data for gurobi python script
y <- seg[sample==standard,c('segment','logR'),with=F]
X <- seg[sample!=standard,c('segment','sample','logR'),with=F]
dat <- merge(y, X, by='segment', all=T)
names(dat) <- c('segment','y','sample','x')
write_tsv(dat,file='/Users/alexgorelick/lab_repos/lpASCN/C161_data_for_gurobi_long.tsv')




mb <- fread(file.path(datadir,'mb_values.tsv'))

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


seg <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_segments.tsv'))
seg <- merge(seg, mb_values[w==0,], by=c('sample'), all.x=T)
seg <- seg[order(sample, segment),]
seg[,logRadj:=logR*m - b]
seg[,global_seg_start:=seg_start/1e6 + global_start_mb]
seg[,global_seg_end:=seg_end/1e6 + global_start_mb]

tmp1 <- seg[Chromosome %in% c(1:22,'X'),c('sample','segment','global_seg_start','global_seg_end','logR'),with=F]
tmp1$type <- 'Original'
tmp2 <- seg[Chromosome %in% c(1:22,'X'),c('sample','segment','global_seg_start','global_seg_end','logRadj'),with=F]
setnames(tmp2,'logRadj','logR')
tmp2$type <- 'Aligned'
plotdat <- rbind(tmp1, tmp2)
plotdat <- merge(plotdat, same_values[w==0,c('segment','same'),with=F], by='segment', all.x=T)
plotdat$type <- factor(plotdat$type, levels=c('Original','Aligned'))
highlight <- plotdat[sample==standard & same==1 & !is.na(same)]
seg$same <- seg$segment %in% highlight$segment
chr <- get_chr_arms()$chr
chr <- chr[chr %in% c(1:22,'X')]
tumor_samples <- unique(seg$sample)
tumor_samples <- tumor_samples[tumor_samples!=standard]
cols <- c('black',rainbow(length(tumor_samples)))
names(cols) <- c(standard, tumor_samples)
p <- ggplot(plotdat) +
    scale_x_continuous(breaks=chr$global_midpoint, labels=chr$chr, expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    geom_vline(xintercept=chr$global_end, linewidth=0.25, color='black') +
    annotate("rect", xmin=highlight$global_seg_start, xmax=highlight$global_seg_end, ymin=min(plotdat$logR)-0.25, ymax=max(plotdat$logR)+0.25, fill='steelblue', alpha=0.03) +
    geom_segment(data=plotdat[sample!=standard], aes(x=global_seg_start, xend=global_seg_end, y=logR, yend=logR, color=sample), linewidth=0.5, alpha=0.25) +
    geom_segment(data=plotdat[sample==standard], aes(x=global_seg_start, xend=global_seg_end, y=logR, yend=logR, color=sample), linewidth=1) +
    scale_color_manual(values=cols, name='Sample') +
    facet_wrap(facets=~type, ncol=1) +
    theme_fit(base_size=10) +
    labs(x='Genomic position', y='logR')
ggsave('aligned_logR.pdf',width=11, height=7)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check the 'same' segments in each sample's adjusted logR, which ones are NOT the same?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ape)
library(ggtree)
library(cowplot)

qc <- seg[Chromosome %in% c(1:22,'X')]
qc <- data.table::dcast(segment ~ sample, value.var='logRadj', data=qc)
aligned_seg <- same_values[same==1,(segment)]
qc <- d2m(qc)
diff_from_standard <- qc - qc[,standard]

## hierarchical clustering of samples based on logRadj difference from standard
dm <- dist(t(diff_from_standard))
tree <- phangorn::upgma(dm)
p1 <- ggtree(tree)
tmp <- as.data.table(p1$data)
tmp <- tmp[isTip==T,]
tmp <- tmp[order(y),]

## get segment order
segment_order <- c(as.integer(unique(pd[Aligned=='No',(segment)])),
                   as.integer(unique(pd[Aligned=='Yes',(segment)])))

pd <- cbind(segment=as.integer(rownames(diff_from_standard)), as.data.table(diff_from_standard))
pd <- data.table::melt(pd, id.var='segment')
names(pd) <- c('segment','sample','diff')
pd$segment <- factor(pd$segment, levels=sort(unique(pd$segment))) #segment_order)
pd[segment %in% aligned_seg, Aligned:='Aligned']
pd[!segment %in% aligned_seg, Aligned:='Not aligned']
pd$Aligned <- factor(pd$Aligned, levels=c('Not aligned','Aligned'))
pd$sample <- factor(pd$sample, levels=tmp$label)

p2 <- ggplot(pd, aes(x=segment, y=sample)) +
    scale_y_discrete(position='right') +
    geom_tile(aes(fill=diff), color=NA) +
    scale_fill_gradient2(low='blue', mid='white', high='red', name='Diff from "standard" logR') +
    facet_grid (.~ Aligned, scales = "free_x", space = "free_x") +
    theme_fit(base_size=12) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position='bottom') +
    labs(x='Segment', y='Sample')

p <- plot_grid(p1, p2, align='h', ncol=2, rel_widths=c(1,5), axis=c('tb'))
ggsave('logR_diff_from_standard_heatmap.pdf',width=11,height=7)


source('func.R')

seg <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_segments.tsv'))
d <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))
standard <- 'PT8-A'
sex <- 'female'

## first try fitting based on dipLogR?
## goal is to get 2q, 11p, 14q, 19p on integer copies
obj <- get_obj_for_sample(standard, seg, d, sex=sex)
fit <- get_fit(obj, dipLogR=-0.46, cores=4)
plot_fit(fit, obj, highlight_seg=aligned_seg)

## not happy with the negative CNs and aligned segs not all close to integers, trying manually
fit <- get_fit(obj, purity=0.75, ploidy=4.4, cores=4)
plot_fit(fit, obj, int_copies=F, highlight_seg=aligned_seg)






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scrap below
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples <- unique(seg$sample)
cores <- 6
this.sex <- 'female'
this.normal <- 'Normal1'

dipLogRs <- seq(-2,2,by=0.01)
cores <- 6
this.sex <- 'female'
this.normal <- 'Normal1'
l <- sweep_sample(tumor_sample,  sex, seg, d, dipLogRs, cores)



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




