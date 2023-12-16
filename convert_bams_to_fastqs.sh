#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem 4000
#SBATCH -o run_bam2fq_%A_%a.out
#SBATCH --array=0-17
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

# This pipeline should be run from FASTQ files to ensure that all preprocessing is standardized. In case you only have access to BAM files and you aren't sure how they were generated, you can convert them back to FASTQs with this script.

module load picard/2.27.5
mkdir -p 00_fastq

## declare arrays
declare -a samples=("C157B1" "C157B2" "C157Di1" "C157Di2" "C157LN4" "C157N1" "C157P10" "C157P1" "C157P2" "C157P3" "C157P4" "C157P6" "C157P9" "C157TD1" "C157TD2a" "C157TD2b" "C157TD3" "C157TD7")

i=${SLURM_ARRAY_TASK_ID}
sample=${samples[$i]}

## convert the bam file back to PE-read FASTQs
BAM="original_bams/${sample}_aligned.bam"
FQ1="00_fastq/${sample}_R1_001.fastq"
FQ2="00_fastq/${sample}_R2_001.fastq"
java -jar $PICARD/picard.jar SamToFastq I=$BAM FASTQ=$FQ1 SECOND_END_FASTQ=$FQ2

## compress the fastq files
gzip $FQ1
gzip $FQ2
