# lpASCN
scripts for generating allele-specific copy number data for multi-region bulk tumor samples with low-pass WGS

## Obtaining het-SNPs from lpWGS data
Azenta's lpWGS analysis pipeline includes genotype imputation: https://web.azenta.com/genotyping-arrays-low-pass-genome-sequencing. These are provided in files named `analysis/[samplename].vcf.gz` (with genome build b37). You can subset these for only high-confidence heterozygous SNPs via:

```
bcftools view -i "%FILTER='PASS'" -g het -O z [samplename].vcf.gz > [samplename]_hetsnps.vcf.gz
```
