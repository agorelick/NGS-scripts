#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem 12000
#SBATCH -o run_filtermutectcalls.out
#SBATCH --mail-user=sebastian_degner@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL


## load modules for bwa/gatk
module load gcc/6.2.0 gatk/4.1.9.0 bcftools/1.13

reference="/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta"
tmp_dir="/n/scratch/users/s/sed362"

## input files
unfiltered_vcf="mutect2/SD8_mutect2.vcf"
f1r2_file="mutect2/SD8_mutect2_f1r2.tar.gz"

## output files
artifact_prior_tables="mutect2/SD8_artifact-prior.tar.gz"
filtered_vcf="mutect2/SD8_filtered.vcf.gz"
passed_vcf="mutect2/SD8_filtered_passed.vcf.gz"



## learn read orientation bias
if test -f $artifact_prior_tables; then
    echo "${artifact_prior_tables} already exists!"
else
    echo "Running LearnReadOrientationModel."
    gatk LearnReadOrientationModel --input $f1r2_file --output $artifact_prior_tables
fi


## filter mutect calls
if test -f $filtered_vcf; then
    echo "${filtered_vcf} already exists!"
else
    echo "Running FilterMutectCalls."
    gatk FilterMutectCalls -R $reference -V $unfiltered_vcf --orientation-bias-artifact-priors $artifact_prior_tables -O $filtered_vcf

    echo "Running IndexFeatureFile."
    gatk IndexFeatureFile -I $filtered_vcf --tmp-dir $tmp_dir
fi


## subsetting passed mutect calls
if test -f $passed_vcf; then
    echo "${passed_vcf} already exists!"
else
    echo "Subsetting for PASS mutations."
    bcftools view -i"FILTER='PASS'" $filtered_vcf | bcftools view -I -O z -o $passed_vcf -

    echo "Running IndexFeatureFile."
    gatk IndexFeatureFile -I $passed_vcf --tmp-dir $tmp_dir
fi

~
