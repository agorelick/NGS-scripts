run_multipcf <-
function(d, sex, build='hg19', penalty=70, tmpdir='.', cleanup=T, output_dir='./multipcf', seed=as.integer(Sys.time())) {

    message('Preparing to run multipcf from ASCAT.')
    t_logR <- data.table::dcast(bin + Chromosome + Position ~ sample, data=d, value.var='logR')
    t_logR[,Normal1:=NULL]
    t_logR <- d2m(t_logR)
    n_logR <- dcast(bin + Chromosome + Position ~ sample, data=d[sample %in% c('Normal1'),], value.var='logR')
    n_logR <- d2m(n_logR)
    t_BAF <- data.table::dcast(bin + Chromosome + Position ~ sample, data=d, value.var='BAF')
    t_BAF[,Normal1:=NULL]
    t_BAF <- d2m(t_BAF)
    n_BAF <- dcast(bin + Chromosome + Position ~ sample, data=d[sample %in% c('Normal1'),], value.var='BAF')
    n_BAF <- d2m(n_BAF)

    if(!dir.exists(tmpdir)) dir.create(tmpdir)
    tumor_logR_file <- file.path(tmpdir,'Tumor_LogR.txt')
    germline_logR_file <- file.path(tmpdir,'Germline_LogR.txt')
    tumor_BAF_file <- file.path(tmpdir,'Tumor_BAF.txt')
    germline_BAF_file <- file.path(tmpdir,'Germline_BAF.txt')
    write_distance_matrix(t_logR,tumor_logR_file)
    write_distance_matrix(n_logR,germline_logR_file)
    write_distance_matrix(t_BAF,tumor_BAF_file)
    write_distance_matrix(n_BAF,germline_BAF_file)

    ascat.bc = ascat.loadData(Tumor_LogR_file = tumor_logR_file,
                              Tumor_BAF_file = tumor_BAF_file,
                              Germline_LogR_file = germline_logR_file,
                              Germline_BAF_file = germline_BAF_file,
                              gender = rep(sex,nrow(t_logR), genomeVersion = build), chrs=c(1:22,'X','Y','MT')) 
    message('Running ascat.asmultipcf() with penalty=',penalty,', seed=',seed,'.')
    ascat.bc = ascat.asmultipcf(ascat.bc, out.dir=NA, seed=seed, penalty=penalty)

    if(cleanup==T) {
        message('Cleaning up after multipcf.')
        trash <- file.remove(c(tumor_logR_file,germline_logR_file,tumor_BAF_file,germline_BAF_file))
    }


    message('Writing extracting segmented BAF and logR data from multipcf for each sample to ', output_dir,'.')
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
    multipcf.extracted$Chromosome <- factor(multipcf.extracted$Chromosome, levels=c(1:22,'X','Y','MT'))
    multipcf.extracted <- multipcf.extracted[order(sample, Chromosome, Position),]
    multipcf.extracted$bin <- as.integer(gsub('bin','',multipcf.extracted$bin)) - 1
    if(!dir.exists(output_dir)) dir.create(output_dir)

    for(s in unique(multipcf.extracted$sample)) {
        message(s)
        write_tsv(multipcf.extracted[sample==s,], file.path(output_dir,paste0(s,'_multipcf.txt')))
    }
}
