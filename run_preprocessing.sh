#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem 16000
#SBATCH -o run_preprocessing_%A_%a.out
#SBATCH --array=0-8   # run a different job for each sample (here for 9 samples, starting at 0)
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

## Notes:
## this script generated analysis-ready bam files based on GATK's Best Practices for Data pre-processing for variant discovery (see: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912).
## run on o2: "sbatch run_preprocessing.sh"
## run this script from a directory containing subdirectory "00_fastq/", which contains paired-end FASTQ files with suffix: "_R1_001.fastq.gz" and "_R2_001.fastq.gz". The adapter sequences in the cutadapt command are for Illumina universal adapter sequences, which is what Azenta uses for lpWGS. Double check this by running FASTQC on your FASTQ files.
## this script is setup for alignment to human reference genome b37 (Broad's version of GRCh37/hg19). This can be changed to b38/GRCh38/hg38, but as of 12/16/2023 Alex's haplotype-aware copy number calling pipeline requires b37.

mkdir -p cutadapt
mkdir -p preprocessing

# use BWAIndex/version0.6.0 (needs to be greater than 0.5.* for bwa>=0.7)
reference='/home/alg2264/data/alex/reference_data/assemblies/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta'
polymorphic_sites='/home/alg2264/data/alex/reference_data/dbSNP/dbSNP_GRCh37/00-common_all.vcf.gz'

## declare vials array (names of tubes submitted to azenta) and samples array (names of samples)
declare -a vials=("1A" "1B" "1C" "1D" "1E" "1F" "1G" "1H" "2A")
declare -a samples=("RLL2" "RML1b" "RLL1" "Normal1" "RUL1c" "Th1" "Kn5" "Kn8" "Kn3")

## get the current sample/read group
i=${SLURM_ARRAY_TASK_ID}
vial=${vials[$i]}
sample=${samples[$i]}
readgroup=$i

## load modules required for cutadapt
module load gcc/9.2.0 python/3.8.12 cutadapt/4.1

if [ ! -f cutadapt/trimmed.${sample}_R1_001.fastq ] || [ ! -f cutadapt/trimmed.${sample}_R2_001.fastq ] ; then
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --cores=20 \
        -o cutadapt/trimmed.${sample}_R1_001.fastq -p cutadapt/trimmed.${sample}_R2_001.fastq \
        00_fastq/${vial}_R1_001.fastq.gz 00_fastq/${vial}_R2_001.fastq.gz
fi


## load modules for bwa/gatk
module load gcc/6.2.0 bwa/0.7.15 gatk/4.1.9.0

if [ ! -f preprocessing/${sample}_raw.sam ]; then
    cmd="bwa mem -M -t 20 -R '@RG\tID:"${i}"\tSM:"${sample}"\tPL:Illumina' $reference cutadapt/trimmed.${sample}_R1_001.fastq cutadapt/trimmed.${sample}_R2_001.fastq > preprocessing/${sample}_raw.sam"
    echo $cmd
    eval $cmd
fi

if [ -f preprocessing/${sample}_raw.sam ] && [ ! -f preprocessing/${sample}_marked_dup.bam ]; then
    cmd="gatk MarkDuplicatesSpark --input preprocessing/${sample}_raw.sam --output preprocessing/${sample}_marked_dup.bam --tmp-dir ~/scratch --reference $reference -M preprocessing/${sample}_marked_dup_metrics.txt --conf 'spark.executor.cores=20'"
    echo $cmd
    eval $cmd
fi

if [ -f preprocessing/${sample}_marked_dup.bam ] && [ ! -f preprocessing/${sample}_recal_data.table ]; then
    cmd="gatk BaseRecalibrator -I preprocessing/${sample}_marked_dup.bam -R $reference --known-sites $polymorphic_sites -O preprocessing/${sample}_recal_data.table --tmp-dir ~/scratch"
    echo $cmd
    eval $cmd
fi

if [ -f preprocessing/${sample}_recal_data.table ] && [ ! -f preprocessing/${sample}_recal.bam ]; then
    cmd="gatk ApplyBQSR -I preprocessing/${sample}_marked_dup.bam -R $reference --bqsr-recal-file preprocessing/${sample}_recal_data.table -O preprocessing/${sample}_recal.bam --tmp-dir ~/scratch"
    echo $cmd
    eval $cmd
fi


