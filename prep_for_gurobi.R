library(here)
library(data.table)

files <- dir('~/Dropbox (Partners HealthCare)/MGH CLINICAL/ALEX/germline_data/lpASCN/C161/processed_data/',full.names=T)
files <- grep('[.]rds',files,value=T)
l <- lapply(files, readRDS)
d <- rbindlist(l)
d$nai <- round(d$na)
d$nbi <- round(d$nb)
d[!is.na(nai) & !is.na(nbi), cn:=nai+nbi]
d[!is.na(nai) & is.na(nbi), cn:=nai]
d[is.na(nai) & is.na(nbi), cn:=-999]
out <- d[,c('sample','segment','dipLogR','cn'),with=F]
names(out) <- c('sample','segment','profile','cn')
write_tsv(out,'C161_data_for_gurobi.tsv')
