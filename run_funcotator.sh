#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem 8000
#SBATCH -o run_funcotator.out
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

module load gcc/6.2.0 gatk/4.1.9.0
mkdir -p somatic_maf

## reference data
reference="/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta"
datasources="/n/data1/hms/genetics/naxerova/lab/alex/reference_data/funcotator_dataSources.v1.7.20200521s"
tmp_dir="/home/alg2264/scratch"

## variable file paths
input_vcf=mutect2/C157_filtered_passed.vcf.gz
overall_maf=somatic_maf/C157_filtered_passed.maf
anno_maf=somatic_maf/C157_filtered_passed_anno.maf
 
## declare sample array (must match order in the VCF file!)
declare -a samples=("C157B1" "C157B2" "C157Di1" "C157Di2" "C157LN4" "C157N1" "C157P10" "C157P1" "C157P2" "C157P3" "C157P4" "C157P6" "C157P9" "C157TD1" "C157TD2a" "C157TD2b" "C157TD3" "C157TD7")



## use funcotator on the passed mutations. However, this will be missing sample-specific data like read-depths. Those are added below.
if test -f $overall_maf; then
    echo "${overall_maf} already exists!"
else
    echo "Indexing the PASS VCF"
    gatk IndexFeatureFile -I $input_vcf

    echo "Running gatk Funcotator."
    gatk Funcotator \
        --variant $input_vcf \
        --reference $reference \
        --ref-version hg19 \
        --data-sources-path $datasources \
        --output $overall_maf \
        --output-file-format MAF \
        --remove-filtered-variants false \
        --tmp-dir $tmp_dir \
        --add-output-vcf-command-line false \
        --transcript-selection-mode CANONICAL \
        --annotation-override Center:naxerovalab \
        --annotation-override Tumor_Sample_Barcode:NA
fi


## Merge in the sample-specific data to the annotated MAF file
module load gcc/9.2.0 bcftools/1.14
tmp_data=${input_vcf}.tmp.tsv   # temp file used for merging
tmp_data_to_merge=${input_vcf}_for_merging.tmp.tsv   # temp file used for merging


if test -f $tmp_data; then
    echo "${tmp_data} already exists!"
else
    echo "Running bcftools query."

    ## create the header for the tmp file
    fields=''
    for sample in "${samples[@]}"
    do
        fields=$fields"\tGT_${sample}\tAD_${sample}\tAF_${sample}\tDP_${sample}\tF1R2_${sample}\tF2R1_${sample}\tSB_${sample}\tPGT_${sample}\tPID_${sample}"
    done
    fields="_CHROM\t_POS\t_REF\t_ALT\tFILTER"${fields}
    echo -e $fields > $tmp_data

    # extract relevant data from VCF which is lost by GATK Funcotator to the tmp file
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%AD\t%AF\t%DP\t%F1R2\t%F2R1\t%SB\t%PGT\t%PID]\n' $input_vcf >> $tmp_data
fi


if test -f $anno_maf; then
    echo "${anno_maf} already exists!"
else
    echo "Joining annotated MAF with data extracted from VCF."

    ## start file with the new header
    header1="$(grep -v '^#' ${overall_maf} | head -1)"
    header2="$(grep -v '^#' ${tmp_data} | head -1)"
    newheader=${header1}" "${header2}
    echo $newheader | awk '{ gsub(" ", "\t") ; print $0 }' > $anno_maf

    ## append joined data
    join -t ':' -1 1 -2 1  <(awk -F '\t' '{print $5$6$11$13 ":" $0;}' $overall_maf | sort -t ':' -k1,1 ) <(awk -F '\t' '{print $1$2$3$4 ":SPLITCOLHERE" $0;}' $tmp_data | sort -t ':' -k1,1 ) | cut -d ':' -f 2- >> $anno_maf

    ## replace temporary string at joining with \t
    sed -i 's/:SPLITCOLHERE/\t/g' $anno_maf

    ## compress the resulting maf
    gzip $anno_maf
fi


