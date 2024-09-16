library(ASCAT)

args = as.character(commandArgs(trailingOnly=TRUE))
patient <- args[1]
tumor_bam <- args[2]
tumor_sample <- args[3]
normal_bam <- args[4]
normal_sample <- args[5]
sex <- args[6]

## output filenames
tumourLogR_file = paste(patient,tumor_sample,normal_sample,"Tumor_LogR.txt",sep='_')
tumourBAF_file = paste(patient,tumor_sample,normal_sample,"Tumor_BAF.txt",sep='_')
normalLogR_file = paste(patient,tumor_sample,normal_sample,"Germline_LogR.txt",sep='_')
normalBAF_file = paste(patient,tumor_sample,normal_sample,"Germline_BAF.txt",sep='_')

if(!file.exists(tumourLogR_file)) {
    ascat.prepareHTS(
                     tumourseqfile = tumor_bam,
                     tumourname = tumor_sample,
                     normalseqfile = normal_bam,
                     normalname = paste0(tumor_sample,'_',normal_sample),
                     allelecounter_exe = "/home/alg2264/.conda/envs/alleleCount/bin/alleleCounter",
                     alleles.prefix = "/home/alg2264/data/alex/reference_data/ascat/G1000_allelesAll_hg19/G1000_alleles_hg19_chr",
                     loci.prefix = "/home/alg2264/data/alex/reference_data/ascat/G1000_lociAll_hg19/G1000_loci_hg19_chr",
                     gender = sex,
                     genomeVersion = "hg19",
                     nthreads = 4,
                     tumourLogR_file = tumourLogR_file,
                     tumourBAF_file = tumourBAF_file,
                     normalLogR_file = normalLogR_file,
                     normalBAF_file = normalBAF_file,
                     BED_file = "/home/alg2264/data/alex/reference_data/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list.bed"
                     )
} else {
    stop('Output already exists!')
}


