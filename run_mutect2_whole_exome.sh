#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-48:00
#SBATCH -p medium
#SBATCH --mem 30000
#SBATCH -o run_mutect2.out
#SBATCH --mail-user=@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

## load required modules
module load gcc/6.2.0 bwa/0.7.15 gatk/4.1.9.0

# use BWAIndex/version0.6.0 (needs to be greater than 0.5.* for bwa>=0.7)
reference='/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta'
intervals_targets='/n/data1/hms/genetics/naxerova/lab/alex/reference_data/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list'
PoN='/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/Mutect2-exome-panel.vcf'
germline_resource='/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.b37/af-only-gnomad.raw.sites.b37_exome.vcf'

# define the output files
output_vcf='SD8_mutect2.vcf'
output_f1r2='SD8_mutect2_f1r2.tar.gz'

# replace this with your directory for large temporary files.
tmp_dir='~/scratch'

# Run Mutect2
if [ -f "$output_f1r2" ]; then
    echo "$output_f1r2 exists. Not regenerating it."
else
    echo "Running Mutect2 ...\n"
    gatk Mutect2 \
        --tmp-dir $tmp_dir \
        --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx30G ' \
        -R $reference \
        -I #/data/naxerova/C157_WES/30-860665034/processed/C157B1_sorted_markdup_withRG_recal.bam \
        -I #/data/naxerova/C157_WES/30-860665034/processed/C157B2_sorted_markdup_withRG_recal.bam \
        -I preprocessing/SD8_NC1-1-1_recal.bam \
        -L $intervals_targets \
        -normal "SD8_NC1-1-1" \
        --germline-resource $germline_resource \
        --panel-of-normals $PoN \
        -O $output_vcf \
        --f1r2-tar-gz $output_f1r2
fi
