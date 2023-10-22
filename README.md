# lpASCN
scripts for generating allele-specific copy number data for multi-region bulk tumor samples with low-pass WGS

## Obtaining het-SNPs from lpWGS data analyzed by Azenta
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

## Filter phased VCF for heterozygous SNPs with high confidence (GP>=0.9)
```
bcftools view -O b -g het C157N1_Normal1_ligated.bcf -i 'FORMAT/GP>=0.90' > C157N1_Normal1_ligated_hetsnps.bcf
```

## Run MEDICC2 for copy-number event based phylogenies

```
medicc2 --bootstrap-nr 200 --bootstrap-method chr-wise --events --plot both MDA11_q=10_Q=20_P=0_1000kbs_segments_for_medicc2.tsv medicc2_bs=200_chr_events
```
