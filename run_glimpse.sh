#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem 20000
#SBATCH -o run_glimpse_%A_%a.out
#SBATCH --array=1-22
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

## Notes:
## This script is run in parallel for each chromosome among 1-22. The reference data 
## assumes you're using b37 as your reference genome. If this is ok, you can use 
## Alex's REF and CHUNKS locations provided. Otherwise, follow the instructions 
## on https://odelaneau.github.io/GLIMPSE/ to create your own reference data.
##
## This script uses the GLIMPSE2 software (https://odelaneau.github.io/GLIMPSE/) 
## to take a *normal* sample's low-coverage WGS bam file, impute heterozygous common SNPs 
## across the entire genome, and then statistically phase them. For 1x WGS coverage, 
## >80% accuracy for common SNPs in >1% of the population is reasonable using 1000Genomes 
## as a reference panel (as I do currently). This seems to work for copy number purposes, 
## but it can be improved using larger reference panels (see the GLIMPSE2 paper/website).  
##
## Make sure you use a bam file with base-score recalibration following GATK-best practices 
## for best results (this is the "*recal.bam" generated by run_preprocessing.sh).

mkdir -p glimpse

chr=${SLURM_ARRAY_TASK_ID}
BAM="preprocessing/C157N1_recal.bam"
OUTPUT="glimpse/C157N1"

REF="/n/data1/hms/genetics/naxerova/lab/alex/reference_data/GLIMPSE/reference_panel/split/1000GP.chr${chr}"
CHUNKS="/n/data1/hms/genetics/naxerova/lab/alex/reference_data/GLIMPSE/chunks/chunks.chr${chr}.txt"

while IFS="" read -r LINE || [ -n "$LINE" ];
do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    CHR=$(echo ${LINE} | cut -d" " -f2)
    REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
    REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
    /home/alg2264/.GLIMPSE2/GLIMPSE2_phase_static --bam-file ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUTPUT}_${CHR}_${REGS}_${REGE}.bcf --threads 20
done < ${CHUNKS}
