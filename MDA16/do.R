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
datadir <- '~/Dropbox (Partners HealthCare)/MGH CLINICAL/ALEX/germline_data/lpASCN/MDA16'
phased_bcf <- file.path(datadir,'original_data/MDA16_Normal1_hetsnps_phased.bcf')
pileup_data <- file.path(datadir,'original_data/MDA16.pileup')
qdnaseq_data <- file.path(datadir,paste0('original_data/MDA16_1000kbp_withXYMT.rds'))
sample_map <- file.path(datadir,'original_data/MDA16_sample_map.txt') ## table with columns: pileup_order, sample_name, proper_sample_name, qdnaseq_name

#pileup <- fread(pileup_data)
#qdnaseq <- readRDS(file.path(datadir,paste0('original_data/MDA16_1000kbp_withXYMT.rds')))

d <- preprocess_data(qdnaseq_data, pileup_data, phased_bcf, sample_map)
write_tsv(d,file=file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))
seg <- get_segments(d, normal_sample='Normal1', sex='female', penalty=20)
write_tsv(seg,file=file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments.tsv'))

## try clustering BAFs to remove subclones from BAF-fit
refined <- refine_segments(d, seg, k.max=10, normal_sample='Normal1', sex='female', penalty=20)
write_tsv(refined$d,file=file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_preprocessed_refined.tsv'))
write_tsv(refined$seg,file=file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments_refined.tsv'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# checking how clustering works to remove subclones from BAF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#source('~/lab_repos/lpASCN/func.R')
#d1 <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))
#seg1 <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_segments.tsv'))
#obj1 <- get_obj_for_sample('PT8-A', seg1, d1, sex='female')
#fit1 <- get_fit(obj1, purity=0.75, ploidy=4.4, cores=4)
#p1 <- plot_fit(fit1, obj1)

#d2 <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_preprocessed_refined.tsv'))
#seg2 <- fread(file.path(datadir,'processed_data/C161_q=10_Q=20_P=0_1000kbs_segments_refined.tsv'))
#obj2 <- get_obj_for_sample('PT8-A', seg2, d2, sex='female')
#fit2 <- get_fit(obj2, purity=0.75, ploidy=4.4, cores=4)
#p2 <- plot_fit(fit2, obj2)

#p <- plot_grid(p1, p2, ncol=2)
#ggsave('PT8-A_fit_refined_side_by_side.pdf',width=16, height=10)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first use auto-fit function for each sample
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('~/lab_repos/lpASCN/func.R')
d <- fread(file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))
seg <- fread(file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments.tsv'))
tumor_samples <- unique(seg$sample)
autofit_dir <- file.path(datadir,'autofits')

autofit <- function(sample, d, seg, sex, cores) {
    message(sample)
    obj <- get_obj_for_sample(sample, seg, d, sex=sex)
    fit <- get_fit(obj, cores=cores)
    p <- plot_fit(fit, obj)
    ggsave(file.path(autofit_dir,paste0(sample,'_autofit.pdf')),width=9,height=8.5)
}
if(!dir.exists(autofit_dir)) dir.create(autofit_dir)
trash <- lapply(tumor_samples, autofit, d, seg, sex='female', cores=4)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run multipcf for multi-sample segmentation
# then extract aligned segments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('~/lab_repos/lpASCN/func.R')
d <- fread(file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))
seg <- fread(file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments.tsv'))
#dat <- d2m(data.table::dcast(segment ~ sample, value.var='logR', data=seg))
#write_tsv(dat,file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments.tsv'))

## choose our 'standard' sample
standard <- 'PT4'
valid_segs <- seg[seg_length >= 1e7 & sample==standard] ## subset 
#valid_segs <- valid_segs[order(seg_length,decreasing=T),]
#valid_segs <- valid_segs[!duplicated(Chromosome),]
seg <- seg[segment %in% valid_segs$segment,]
seg <- seg[!sample %in% c('LN1','Liv1'),]

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
quantile(abs(tmp.d$diff_from_seg),0.20)  # 0.02864416

## export data for gurobi python script
y <- seg[sample==standard,c('segment','logR'),with=F]
X <- seg[sample!=standard,c('segment','sample','logR'),with=F]
dat <- merge(y, X, by='segment', all=T)
names(dat) <- c('segment','y','sample','x')
write_tsv(dat[!sample %in% c('LN1','Liv1'),],file=file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_data_for_gurobi_long_no_Liv1_LN1.tsv'))

#write_tsv(dat,file=file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_data_for_gurobi_long.tsv'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# continue after LP optimization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mb <- fread(file.path(datadir,'processed_data/mb_values.tsv'))

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


seg <- fread(file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments.tsv'))
seg <- merge(seg, mb_values[w==0,], by=c('sample'), all.x=T)
seg <- seg[order(sample, segment),]
seg[,logRadj:=logR*m - b]
seg[,global_seg_start:=seg_start/1e6 + global_start_mb]
seg[,global_seg_end:=seg_end/1e6 + global_start_mb]
seg <- seg[!is.na(w),]

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
    annotate("rect", xmin=highlight$global_seg_start, xmax=highlight$global_seg_end, ymin=min(plotdat$logR)-0.25, ymax=max(plotdat$logR)+0.25, fill='steelblue', alpha=0.05) +
    geom_segment(data=plotdat[sample!=standard], aes(x=global_seg_start, xend=global_seg_end, y=logR, yend=logR, color=sample), linewidth=0.5, alpha=0.25) +
    geom_segment(data=plotdat[sample==standard], aes(x=global_seg_start, xend=global_seg_end, y=logR, yend=logR, color=sample), linewidth=1) +
    scale_color_manual(values=cols, name='Sample') +
    facet_wrap(facets=~type, ncol=1) +
    theme_fit(base_size=10) +
    labs(x='Genomic position', y='logR')
ggsave(file.path(datadir,'figures/MDA16_q=10_Q=20_P=0_1000kbs_aligned_logR.pdf'),width=11, height=7)

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

pd <- cbind(segment=as.integer(rownames(diff_from_standard)), as.data.table(diff_from_standard))
pd[segment %in% aligned_seg, Aligned:='Aligned']
pd[!segment %in% aligned_seg, Aligned:='Not aligned']
pd$Aligned <- factor(pd$Aligned, levels=c('Not aligned','Aligned'))

pd <- data.table::melt(pd, id.var=c('segment','Aligned'))
names(pd) <- c('segment','Aligned','sample','diff')
pd$segment <- factor(pd$segment, levels=sort(unique(pd$segment)))
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
ggsave(file.path(datadir,'figures/MDA16_q=10_Q=20_P=0_1000kbs_logR_diff_from_standard_heatmap.pdf'),width=11, height=7)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# roughly fit the standard, then automatically fit 
# it and the other samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('~/lab_repos/lpASCN/func.R')

seg <- fread(file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments.tsv'))
d <- fread(file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_preprocessed.tsv'))
standard <- 'PT4'
sex <- 'female'
samples <- sort(unique(seg$sample))
qc <- seg[segment %in% aligned_seg & sample==standard,]

fitdir <- file.path(datadir,'fits')
if(!dir.exists(fitdir)) dir.create(fitdir)

## fit PT4 (standard)
## goal is to get aligned segments on integer copies
obj1 <- get_obj_for_sample('PT4', seg, d, sex=sex)
fit1 <- get_fit(obj1, purity=0.65, ploidy=4.2, cores=4)
plot_fit(fit1, obj1, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_PT4_standard.pdf'),width=8, height=10.5)


## for each sample, our goal is for these segments to be roughly on integers
standard_ascn <- get_na_nb(p=fit1$pu, t=fit1$pl, x=seg[sample==standard & segment %in% aligned_seg])

x <- seg[segment %in% aligned_seg,]
x <- merge(x, standard_ascn[,c('segment','na','nb'),with=F], by='segment', all.x=T)
x$na <- round(x$na)
x$nb <- round(x$nb)
setnames(x,c('na','nb'),c('na_target','nb_target'))
x[,tcn_target:=na_target+nb_target]

get_fit_for_sample <- function(this.sample, x, use.tcn=F) {
    message(this.sample)
    tmp <- x[sample==this.sample]
    combos <- as.data.table(expand.grid(pl=seq(1,6,by=0.01), pu=seq(0.05,0.95,by=0.01)))
    combos$sample <- this.sample
    combos <- merge(combos, tmp, by='sample', all=T, allow.cartesian=T)
    combo_ascn <- get_na_nb(combos$pu, combos$pl, combos)
    if(use.tcn==T) {
        combo_ascn[,tcn:=na+nb]
        combo_ascn[,dist_from_int:=abs(tcn - tcn_target)]
    } else {
        combo_ascn[,dist_from_int:=abs(na - na_target) + abs(nb - nb_target)]
    }
    combo_ascn[is.infinite(dist_from_int), dist_from_int:=NA]
    test_combo <- function(combos) {
        mean_dist_from_int_across_segments <- median(combos$dist_from_int,na.rm=T) 
        list(mean_dist_from_int_across_segments=mean_dist_from_int_across_segments)
    }
    tmp <- combo_ascn[,test_combo(.SD), by=c('pu','pl')]
    tmp <- tmp[order(mean_dist_from_int_across_segments,decreasing=F),]
    tmp$sample <- this.sample
    tmp
}
l <- lapply(samples, get_fit_for_sample, x)
result <- rbindlist(l)
result[!duplicated(sample),]


obj1 <- get_obj_for_sample('PT4', seg, d, sex=sex)
fit1 <- get_fit(obj1, purity=0.63, ploidy=4.21, cores=4)
plot_fit(fit1, obj1, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_PT4.pdf'),width=8, height=10.5)

obj2 <- get_obj_for_sample('PT1', seg, d, sex=sex)
fit2 <- get_fit(obj2, purity=0.56, ploidy=4.37, cores=4)
plot_fit(fit2, obj2, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_PT1.pdf'),width=8, height=10.5)

obj3 <- get_obj_for_sample('PT2', seg, d, sex=sex)
fit3 <- get_fit(obj3, purity=0.53, ploidy=4.39, cores=4)
plot_fit(fit3, obj3, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_PT2.pdf'),width=8, height=10.5)

obj4 <- get_obj_for_sample('PT3', seg, d, sex=sex)
fit4 <- get_fit(obj4, purity=0.52, ploidy=4.25, cores=4)
plot_fit(fit4, obj4, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_PT3.pdf'),width=8, height=10.5)

obj5 <- get_obj_for_sample('Lun1a', seg, d, sex=sex)
fit5 <- get_fit(obj5, purity=0.51, ploidy=4.27, cores=4)
plot_fit(fit5, obj5, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_Lun1a.pdf'),width=8, height=10.5)

obj6 <- get_obj_for_sample('Lun1b', seg, d, sex=sex)
fit6 <- get_fit(obj6, purity=0.51, ploidy=4.33, cores=4)
plot_fit(fit6, obj6, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_Lun1b.pdf'),width=8, height=10.5)

obj7 <- get_obj_for_sample('Liv2', seg, d, sex=sex)
fit7 <- get_fit(obj7, purity=0.54, ploidy=4.72, cores=4)
plot_fit(fit7, obj7, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_Liv2.pdf'),width=8, height=10.5)

## low purity
obj8 <- get_obj_for_sample('Liv1', seg, d, sex=sex)
fit8 <- get_fit(obj8, purity=0.05, ploidy=3.92, cores=4)
plot_fit(fit8, obj8, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_Liv1.pdf'),width=8, height=10.5)

## low purity
obj9 <- get_obj_for_sample('LN1', seg, d, sex=sex)
fit9 <- get_fit(obj9, purity=0.79, ploidy=4.92, cores=4)
plot_fit(fit9, obj9, int_copies=F, highlight_seg=aligned_seg)
ggsave(file.path(datadir,'fits/MDA16_LN1.pdf'),width=9, height=10.5)


fits_file <- file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_fits.tsv')
fits <- rbind(fit1, fit2, fit3, fit4, fit5, fit6, fit7) #, fit8, fit9)  # excluding low purity samples
write_tsv(fits, fits_file)
         


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate data for MEDICC2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fits <- fread(fits_file)
get_ascn_for_sample <- function(i, fits, seg) {
    this.sample <- fits$sample[i]
    message(this.sample)
    pu <- fits$pu[i]
    pl <- fits$pl[i]
    sample_seg <- seg[sample==this.sample]
    ascn <- get_na_nb(p=pu, t=pl, x=sample_seg)
    ascn
}
l <- lapply(1:nrow(fits), get_ascn_for_sample, fits, seg)
results <- rbindlist(l)
write_tsv(results,file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments_fit.tsv'))



medicc <- results[Chromosome %in% 1:22,]
medicc[is.na(nb)]
medicc$cn_a <- round(medicc$na)
medicc$cn_b <- round(medicc$nb)
medicc <- medicc[,c('sample','Chromosome','seg_start','seg_end','cn_a','cn_b'),with=F]
names(medicc) <- c('sample_id', 'chrom', 'start', 'end', 'cn_a', 'cn_b')
medicc[,start:=start/1e6]
medicc[,end:=end/1e6]
medicc$start <- round(medicc$start)
medicc$end <- round(medicc$end)
write_tsv(medicc,file.path(datadir,'processed_data/MDA16_q=10_Q=20_P=0_1000kbs_segments_for_medicc2.tsv'))



results[!is.na(na) & is.na(nb), tcn:=na]
results[!is.na(na) & !is.na(nb), tcn:=na+nb]
results[Chromosome=='MT']
mat <- data.table::dcast(Chromosome + segment ~ sample, value.var='tcn', data=results)
mat[Chromosome %in% 1:22, Normal:=2]
mat[Chromosome %in% c('X','Y'), Normal:=1]
mat <- mat[Chromosome %in% c(1:22,'X','Y')]
mat[,Chromosome:=NULL]
mat <- d2m(mat)
dm <- dist(t(mat), method='manhattan')
tree <- nj(dm)
tree <- phytools::reroot(tree, node=which(tree$tip.label=='Normal'))




