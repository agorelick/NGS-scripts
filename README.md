# Instructions for low-coverage copy number pipeline on o2


## Clone this repo onto o2

```
# Clone this repo on o2
git clone git@github.com:agorelick/lpASCN.git
```

## (Optional) Convert BAM files back to FASTQ files if necessary.

This pipeline should be run from FASTQ files in order to ensure that all preprocessing is standardized. In case you only have access to BAM files (and you aren't sure exactly how they were generated), you can convert them back to FASTQs with this script. Note: this assumes paired-end reads.

First, copy convert_bams_to_fastqs.sh from this repo to a location on o2 where you can save large data files. Then create a symbolic link to a directory with the bam files you will use: 

```
# create sym-link to the bam directory
ln -s [PATH TO DIRECTORY WITH BAM FILES] original_bams
```

Then, modify the convert_bams_to_fastqs.sh script to include the names of the bam files in the "samples" array. You may also need to change the bam file name in the `BAM="original_bams/${sample}_aligned.bam"` to match how your pre-existing bam files are named.

Finally, run the convert_bams_to_fastqs.sh script to regenerate fastq.gz files the read pairs for each sample. Job will be run very quickly, ~ 5 min per sample.

```
sbatch convert_bams_to_fastqs.sh
```



## Create analysis-ready bam files

This step uses the run_processing.sh slurm script to generate analysis-ready bam files based on GATK's Best Practices for Data pre-processing for variant discovery (see: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912). 

Copy the run_processing.sh script to a location on *o2* containing the subdirectory "00_fastq/" which contains paired-end FASTQ files with suffixes: "_R1_001.fastq.gz" and "_R2_001.fastq.gz" (this is the default provided by Azenta). 

Next, modify the script to use your sample names and vial names (which refer to the name from the fastq file). Also edit the number of jobs in the header `#SBATCH --array=0-8` for your number of samples. These are indexed at 0, so this example is for 9 samples to be run parallel.

Note: the adapter sequences in the cutadapt section are for Illumina universal adapter sequences, which is what Azenta uses for lpWGS. You can check this by running FASTQC on your FASTQ files, and change the adapter sequences if necessary.

Also note: This script is setup for alignment to human reference genome b37 (Broad's version of GRCh37/hg19). This can be changed to b38/GRCh38/hg38, but as of 12/16/2023 Alex's haplotype-aware copy number calling pipeline requires b37.

Total run-time per sample is around 6 hours.

Run the slurm job on o2 via the command `sbatch run_preprocessing.sh`.



## Obtained haplotype-phased heterozygous germline SNPs from a normal sample

### Download pre-built executable files for GLIMPSE2 onto o2
```
mkdir ~/.GLIMPSE2; cd ~/.GLIMPSE2
wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_chunk_static
wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_concordance_static
wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_ligate_static
wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_phase_static
wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_split_reference_static
```




### Obtaining het-SNPs from lpWGS data analyzed by Azenta
Azenta's lpWGS analysis pipeline includes genotype imputation: https://web.azenta.com/genotyping-arrays-low-pass-genome-sequencing. These are provided in files named `analysis/[samplename].vcf.gz` (with genome build b37). You can subset these for only high-confidence heterozygous SNPs via the following bcftools command. Note that we need `-c 0` to be added to make sure that the output has field "AC" required by SHAPEIT4 for phasing.

```
bcftools view -i "%FILTER='PASS'" -g het -O z  -c 0 [samplename].vcf.gz > [samplename]_hetsnps.vcf.gz
bcftools index [samplename]_hetsnps.vcf.gz
```

## Do phasing

## Merge phased BCF files for each chromosome

```
files=''
for chr in {1..22}; do files=$files" "[samplename]_hetsnps_phased_chr${chr}.bcf; done
bcftools concat -o [samplename]_hetsnps_phased.bcf $files
bcftools index [samplename]_hetsnps_phased.bcf
```

## Filter phased VCF for heterozygous SNPs with high confidence (GP>=0.9) and remove FORMAT/DS,GP fields from GLIMPSE2
```
bcftools view -O b -g het C157N1_Normal1_ligated.bcf -i 'FORMAT/GP>=0.90' | bcftools annotate -x FORMAT/DS,FORMAT/GP -O b > C157N1_Normal1_ligated_cleaned.bcf
```

## Filter VCF for mutations with FILTER=={'PASS' or '.'} (after running GATK FilterMutectCalls) 
```
bcftools view -i "%FILTER='PASS' | %FILTER='.'" MDA11_filtered.vcf.gz | bcftools view -I - -O z -o MDA11_filtered_passed.vcf.gz
```

## Run MEDICC2 for copy-number event based phylogenies

```
medicc2 --bootstrap-nr 200 --bootstrap-method chr-wise --events --plot both MDA11_q=10_Q=20_P=0_1000kbs_segments_for_medicc2.tsv medicc2_bs=200_chr_events
```
