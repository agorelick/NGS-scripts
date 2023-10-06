# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(QDNAseq)
library(Rsamtools)
library(data.table)
library(Rcpp)
library(Biobase)
library(pctGCdata)
library(here)
library(ASCAT)

source(here('func.R'))

## which sample to run
args = commandArgs(trailingOnly=TRUE)
sample_index <- as.integer(args[1])
cores <- as.integer(args[2])

## load data
d <- fread('/home/alg2264/data/alex/lpASCN/C161/C161_q=10_Q=20_P=0_1000kbs_preprocessed.tsv')
seg <- fread('/home/alg2264/data/alex/lpASCN/C161/C161_q=10_Q=20_P=0_1000kbs_segments.tsv')
samples <- sort(unique(as.character(seg$sample)))
this_sample <- samples[sample_index]
output_file <- paste0('/home/alg2264/data/alex/lpASCN/C161/C161_q=10_Q=20_P=0_1000kbs_dipLogRsweep_',this_sample,'.rds')

## define dipLogR values to sweep through
dipLogRs <- seq(-2,2,by=0.01)
l <- sweep_sample(this_sample, sex='female', seg, d, dipLogRs, cores=cores)

## save result
message('Saving the result to: ',output_file)
saveRDS(l,file=output_file)


