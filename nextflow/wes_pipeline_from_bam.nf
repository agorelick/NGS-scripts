if (!params.csv) {
    error "Missing required parameter: --csv"
}

// Load the CSV into a list of maps (1 per row)
// params.csv is defined from the nextflow run command via --csv <CSV_FILE>
// CSV columns: sample, orig_bam, orig_bai, order, patient, normal_sample, sex, nextflow_dir, output_dir
def df = file(params.csv)
    .text
    .split('\n')
    .collect { it.trim() }
    .findAll { it } // remove empty lines

def headers = df[0].split(',') // column names
def rows = df[1..-1].collect { line ->
    def fields = line.split(',')
    [ : ].withDefault { null }.tap { map ->
        headers.eachWithIndex { h, i ->
            if (i < fields.size()) {
                map[h] = fields[i]
            } else {
                map[h] = null
            }
        }
    }
}

println "Parsed headers: ${headers}"
println "First row: ${rows[0]}"

// Create tuple channel: (sample, input_bam_dir, order)
// input_bam_dir is the directory containing <sample>.bam and <sample>.bam.bai
def sample_ch = Channel.from( rows.collect { [it.sample, it.orig_bam, it.orig_bai, it.order] } )

// Set scalar params using first row with non-null value
def first_full_row = rows.find { it.patient }
params.patient         = first_full_row.patient
params.normal_sample   = first_full_row.normal_sample
params.sex             = first_full_row.sex
params.nextflow_dir    = first_full_row.nextflow_dir
params.output_dir      = first_full_row.output_dir

// dynamically generated parameters
params.bam_dir        = "${params.output_dir}/${params.patient}/bams"
params.mutect_dir     = "${params.output_dir}/${params.patient}/mutect"
params.maf_dir        = "${params.output_dir}/${params.patient}/mafs"
params.mosdepth_dir   = "${params.output_dir}/${params.patient}/mosdepth"
params.mtbam_dir      = "${params.output_dir}/${params.patient}/mtbams"
params.haplocheck_dir = "${params.output_dir}/${params.patient}/haplocheck"
params.ascat_dir      = "${params.output_dir}/${params.patient}/ascat"
params.fastqc_dir     = "${params.output_dir}/${params.patient}/fastqc"
params.tmp_dir        = "${params.nextflow_dir}/${params.patient}/tmp_files"

/*
 * Additional pipeline parameters (use for all WES data)
 */
// Genome reference files
params.build    = "hg38"
params.mt_label = "chrM"
params.ref_fasta = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
params.ref_amb   = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.amb"
params.ref_ann   = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.ann"
params.ref_bwt   = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.bwt"
params.ref_fai   = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.fai"
params.ref_pac   = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.pac"
params.ref_sa    = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.sa"
params.ref_dict  = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.dict"

// additional reference files with index files
params.polymorphic_sites     = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz"
params.polymorphic_sites_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz.tbi"
params.germline_resource     = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz"
params.germline_resource_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz.tbi"
params.panel_of_normals      = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf"
params.panel_of_normals_idx  = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf.idx"
params.targets_bed           = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38.bed"
params.genome_chunks         = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38_50Mbchunks.csv"

// ascat/CNAlign params
params.allelecounter_exe  = "/home/alg2264/miniconda3/envs/CNalign/bin/alleleCounter"
params.alleles_prefix     = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_allelesAll_hg38/G1000_alleles_hg38_chr"
params.loci_prefix        = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_lociAll_hg38/G1000_loci_GRCh38_chr"
params.gccontentfile      = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/GC_G1000_hg38.txt"
params.replictimingfile   = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/RT_G1000_hg38.txt"

println "params.nextflow_dir: ${params.nextflow_dir}"
println "params.output_dir: ${params.output_dir}"

process MAKE_DIRS {
    tag "mkdirs"
    executor 'local'

    input:
    val bam_dir
    val maf_dir
    val mosdepth_dir
    val mtbam_dir
    val haplocheck_dir
    val ascat_dir
    val tmp_dir
    val mutect_dir
    val fastqc_dir

    output:
    path "mkdir_done.txt"

    """
    mkdir -p ${bam_dir}
    mkdir -p ${maf_dir}
    mkdir -p ${mosdepth_dir}
    mkdir -p ${mtbam_dir}
    mkdir -p ${haplocheck_dir}
    mkdir -p ${ascat_dir}
    mkdir -p ${tmp_dir}
    mkdir -p ${mutect_dir}
    mkdir -p ${fastqc_dir}
    touch mkdir_done.txt
    """
}

/*
 * Convert input BAM to paired-end FASTQ files.
 * The input BAM is expected to be coordinate-sorted (as delivered).
 * We queryname-sort first (required by samtools fastq for correct pairing),
 * then convert to R1/R2 FASTQ.
 */
process BAM_TO_FASTQ {
    tag "$sample"
    cpus 8
    memory '32GB'
    time '4h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), path(orig_bam), path(orig_bai), val(sample_order)
    val make_dirs_ch

    output:
    tuple val(sample), val(sample_order), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz")

    script:
    """
    module load gcc/14.2.0 samtools/1.21

    # queryname-sort the input BAM (required for correct R1/R2 pairing by samtools fastq)
    samtools sort -n -@ 8 -m 3G \
        -o ${sample}_qnsorted.bam \
        ${orig_bam}

    # convert to paired FASTQ; singletons/supplementary reads go to /dev/null
    samtools fastq -@ 8 \
        -1 ${sample}_R1.fastq.gz \
        -2 ${sample}_R2.fastq.gz \
        -0 /dev/null \
        -s /dev/null \
        -n \
        ${sample}_qnsorted.bam

    rm ${sample}_qnsorted.bam
    """
}

/*
 * Run FastQC on raw FASTQs (pre-trimming QC) and trim Illumina Universal Adapters.
 * Combined into a single process so the raw FASTQs can be deleted after both
 * steps are complete — they are consumed by nothing else downstream.
 */
process FASTQC_AND_TRIM {
    tag "$sample"
    cpus 8
    memory '32GB'
    time '4h'
    executor 'slurm'
    queue 'short'
    publishDir params.fastqc_dir, mode: 'copy', pattern: '*_fastqc.{html,zip}'

    input:
    tuple val(sample), val(sample_order), path(fq1), path(fq2)

    output:
    tuple val(sample), val(sample_order), path("${sample}_trimmed_R1.fastq.gz"), path("${sample}_trimmed_R2.fastq.gz"), emit: trimmed
    tuple path("${sample}_R1_fastqc.html"), path("${sample}_R1_fastqc.zip"),
          path("${sample}_R2_fastqc.html"), path("${sample}_R2_fastqc.zip"), emit: fastqc_reports

    script:
    """
    # FastQC on raw (pre-trimming) FASTQs
    conda run -n fastqc_0.11.5 fastqc --threads 2 --outdir . ${fq1} ${fq2}

    # Trim Illumina Universal Adapters
    conda run -n cutadapt cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --minimum-length 20 --cores=8 \
        -o ${sample}_trimmed_R1.fastq.gz \
        -p ${sample}_trimmed_R2.fastq.gz \
        ${fq1} ${fq2}

    # raw FASTQs are no longer needed — both consumers (FastQC and cutadapt) are done
    rm ${fq1} ${fq2}
    """
}

/*
 * Run BWA-MEM
 * Deletes trimmed FASTQs after alignment — they are the sole consumer.
 * The rm targets the staged symlink names, which frees the work-dir copy
 * created by TRIM_ADAPTERS.
 */
process BWA_MEM {
    tag "$sample"
    cpus 8
    memory '32GB'
    time '8h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), val(sample_order), path(trimmed_fq1), path(trimmed_fq2)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_raw.sam")

    script:
    """
    module load gcc/14.2.0 bwa/0.7.18

    bwa mem -M -t 8 -R '@RG\\tID:${sample_order}\\tSM:${sample}\\tPL:Illumina' \
        ${ref_fasta} ${trimmed_fq1} ${trimmed_fq2} > ${sample}_raw.sam

    # trimmed FASTQs are no longer needed after alignment
    rm ${trimmed_fq1} ${trimmed_fq2}
    """
}

/*
 * Run samtools sort
 * Deletes the raw SAM after sorting — SAMTOOLS_SORT is the sole consumer.
 */
process SAMTOOLS_SORT {
    tag "$sample"
    cpus 8
    memory '72GB'
    time '3h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), path(raw_sam)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_sorted.bam"), path("${sample}_sorted.bam.bai")

    script:
    """
    module load gcc/14.2.0 samtools/1.21

    samtools sort -@ 8 -m 3G -o ${sample}_sorted.bam ${raw_sam}
    samtools index -@ 8 ${sample}_sorted.bam

    # raw SAM is no longer needed after sorting
    rm ${raw_sam}
    """
}

/*
 * Run GATK MarkDuplicates
 * Deletes the sorted BAM/BAI after duplicate marking — GATK_MARKDUP is the
 * sole consumer. The publishDir copies realigned.bam to the output directory
 * before Nextflow stages-out, so deleting sorted.bam here is safe.
 */
process GATK_MARKDUP {
    tag "$sample"
    cpus 8
    memory '64GB'
    time '12h'
    executor 'slurm'
    queue 'short'
    publishDir params.bam_dir, mode: 'copy'

    input:
    tuple val(sample), path(sorted_bam), path(sorted_bai)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_realigned.bam"), path("${sample}_realigned.bam.bai"), path("${sample}_realigned_marked_dup_metrics.txt")

    script:
    """
    conda run -n gatk_4.6.1.0 gatk --java-options "-Xmx56g" MarkDuplicatesSpark \
        -I ${sorted_bam} \
        -O ${sample}_realigned.bam \
        -M ${sample}_realigned_marked_dup_metrics.txt \
        -R ${ref_fasta} \
        --create-output-bam-index true \
        --spark-master "local[8]"

    # sorted BAM/BAI are no longer needed after duplicate marking
    rm ${sorted_bam} ${sorted_bai}
    """
}

/*
 * Run GATK Mutect2 for multi-sample tumor/normal variant calling
 */
process GATK_MUTECT2 {
    tag "$region"
    cpus 18
    memory '32GB'
    time '24h'
    executor 'slurm'
    queue 'medium'
    publishDir params.mutect_dir, mode: 'copy'

    input:
    tuple val(chr), val(start), val(end), val(region)
    path bed_file
    path all_bams
    path all_bam_indices
    val patient
    val normal_sample
    path output_dir
    tuple path(polymorphic_sites), path(polymorphic_sites_tbi)
    tuple path(germline_resource), path(germline_resource_tbi)
    tuple path(panel_of_normals), path(panel_of_normals_idx)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(region), path("regions_${region}.bed"), path("${patient}_${region}.vcf.gz"), path("${patient}_${region}.vcf.gz.tbi"), path("${patient}_${region}.vcf.gz.stats"), path("${patient}_${region}.f1r2.tar.gz")

    script:
    def bams_line = all_bams.collect { bam -> "-I ${bam}" }.join(' ')
    """
    module load bedtools/2.31.0

    echo -e "${chr}\t${start}\t${end}" | bedtools intersect -a ${bed_file} -b - > regions_${region}.bed

    conda run -n gatk_4.6.1.0 gatk Mutect2 -R $ref_fasta \
        $bams_line \
        -normal $normal_sample \
        -L regions_${region}.bed \
        --f1r2-tar-gz ${patient}_${region}.f1r2.tar.gz \
        --native-pair-hmm-threads 16 \
        -O ${patient}_${region}.vcf.gz \
        --germline-resource $germline_resource \
        --panel-of-normals $panel_of_normals \
    """
}

/*
 * Merge output from mutect2 (VCFs, f1r2-files, stats-files), across genomic chunks
 * Deletes per-region VCFs, stats, and f1r2 tarballs after merging.
 */
process MERGE_REGIONS {
    tag "$patient"
    cpus 8
    memory '32GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    input:
    val patient
    path all_vcf
    path all_vcf_tbi
    path all_stats
    path all_f1r2
    publishDir params.mutect_dir, mode: 'copy'

    output:
    tuple path("${patient}_raw.vcf.gz"), path("${patient}_raw.vcf.gz.tbi"), path("${patient}_raw.vcf.gz.stats"), path("${patient}_raw.artifact-prior.tar.gz")

    script:
    def vcf_line   = all_vcf.collect   { vcf       -> "${vcf}"            }.join(' ')
    def stats_line = all_stats.collect { statsfile  -> "--stats ${statsfile}" }.join(' ')
    def f1r2_line  = all_f1r2.collect  { f1r2file   -> "-I ${f1r2file}"   }.join(' ')
    """
    module load bcftools/1.21

    bcftools concat ${vcf_line} -a > ${patient}_raw_unsorted.vcf
    bcftools sort ${patient}_raw_unsorted.vcf -O z -o ${patient}_raw.vcf.gz
    rm ${patient}_raw_unsorted.vcf

    conda run -n gatk_4.6.1.0 gatk IndexFeatureFile -I ${patient}_raw.vcf.gz
    conda run -n gatk_4.6.1.0 gatk MergeMutectStats ${stats_line} --output ${patient}_raw.vcf.gz.stats
    conda run -n gatk_4.6.1.0 gatk LearnReadOrientationModel ${f1r2_line} --output ${patient}_raw.artifact-prior.tar.gz

    # per-region intermediates are no longer needed after merging
    rm -f ${all_vcf.join(' ')} ${all_vcf_tbi.join(' ')} ${all_stats.join(' ')} ${all_f1r2.join(' ')}
    """
}

/*
 * Filter Mutect calls
 */
process FILTER_MUTECT_CALLS {
    tag "$patient"
    cpus 8
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'
    publishDir params.mutect_dir, mode: 'copy'

    input:
    val patient
    tuple path(raw_vcf), path(raw_vcf_tbi), path(raw_vcf_stats), path(raw_artifact_priors)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple path("${patient}_unfiltered_norm.vcf.gz"), path("${patient}_unfiltered_norm.vcf.gz.tbi"), path("${patient}_filtered.vcf.gz"), path("${patient}_filtered.vcf.gz.tbi")

    script:
    """
    module load bcftools/1.21

    conda run -n gatk_4.6.1.0 gatk FilterMutectCalls \
        -R $ref_fasta -V $raw_vcf \
        --orientation-bias-artifact-priors $raw_artifact_priors \
        -O ${patient}_unfiltered.vcf.gz

    conda run -n gatk_4.6.1.0 gatk IndexFeatureFile -I ${patient}_unfiltered.vcf.gz

    bcftools norm --multiallelics -both --fasta-ref $ref_fasta \
        ${patient}_unfiltered.vcf.gz \
        | bcftools view -I -O z -o ${patient}_unfiltered_norm.vcf.gz -

    conda run -n gatk_4.6.1.0 gatk IndexFeatureFile -I ${patient}_unfiltered_norm.vcf.gz

    bcftools view -i "FILTER='PASS'" ${patient}_unfiltered_norm.vcf.gz \
        | bcftools view -I -O z -o ${patient}_filtered.vcf.gz -

    conda run -n gatk_4.6.1.0 gatk IndexFeatureFile -I ${patient}_filtered.vcf.gz

    # intermediate unfiltered VCF (pre-normalization) no longer needed
    rm -f ${patient}_unfiltered.vcf.gz ${patient}_unfiltered.vcf.gz.tbi
    """
}

/*
 * Run VCF2MAF on filtered VCF file for each sample
 */
process VCF2MAF {
    tag "$sample"
    cpus 1
    memory '16GB'
    time '1h'
    executor 'slurm'
    queue 'short'
    publishDir params.maf_dir, mode: 'copy'

    input:
    tuple val(sample), val(bam_dir), val(sample_order)
    path filtered_vcf
    path filtered_vcf_tbi

    output:
    path "${sample}.maf"

    script:
    """
    module load bcftools/1.21

    bcftools view $filtered_vcf -s $sample > ${sample}.vcf

    conda run -n vep perl /home/alg2264/repos/vcf2maf/vcf2maf.pl \
        --input-vcf ${sample}.vcf \
        --output-maf ${sample}.maf \
        --tumor-id ${sample} \
        --remap-chain /home/alg2264/repos/vcf2maf/data/hg38_to_GRCh38.chain

    # per-sample VCF slice is no longer needed after MAF conversion
    rm ${sample}.vcf
    """
}

/*
 * Run mosdepth to get bed file with x-fold coverage in targeted regions
 */
process MOSDEPTH {
    tag "$sample"
    cpus 2
    memory '8GB'
    time '2h'
    executor 'slurm'
    queue 'short'
    publishDir params.mosdepth_dir, mode: 'copy'

    input:
    tuple val(sample), path(sample_bam), path(sample_bai)
    path bed_file

    output:
    path "${sample}.mosdepth.summary.txt"
    path "${sample}.mosdepth.global.dist.txt"
    path "${sample}.regions.bed.gz.csi"
    path "${sample}.regions.bed.gz"
    path "${sample}.mosdepth.region.dist.txt"

    script:
    """
    conda run -n mosdepth mosdepth -n -t 2 -m --by ${bed_file} ${sample} ${sample_bam}
    """
}

/*
 * Slice bam for chrM/MT for haplocheck
 */
process SLICE_MTDNA {
    tag "$sample"
    cpus 1
    memory '4GB'
    time '5m'
    executor 'slurm'
    queue 'short'
    publishDir params.mtbam_dir, mode: 'copy'

    input:
    tuple val(sample), path(sample_bam), path(sample_bai)
    val mt_label

    output:
    tuple val(sample), path("${sample}_mt.bam"), path("${sample}_mt.bam.bai")

    script:
    """
    module load gcc/14.2.0 samtools/1.21
    samtools view -hb ${sample}_realigned.bam -F 3840 -q 20 ${mt_label} > ${sample}_mt.bam
    samtools index ${sample}_mt.bam
    """
}

/*
 * Prepare CNA data (allele counts for ASCAT/CNAlign)
 */
process PREP_CNA_DATA {
    tag "$tumor_sample"
    cpus 2
    memory '16GB'
    time '2h'
    executor 'slurm'
    queue 'short'
    publishDir params.ascat_dir, mode: 'copy'

    input:
    tuple val(tumor_sample), path(tumor_bam), path(tumor_index), val(normal_sample), path(normal_bam), path(normal_index)
    val patient
    val sex
    val build
    path bed_file
    path allelecounter_exe
    val alleles_prefix
    val loci_prefix

    output:
    tuple(
        path("${patient}_${tumor_sample}_${normal_sample}_Germline_BAF_rawBAF.txt"),
        path("${patient}_${tumor_sample}_${normal_sample}_Germline_BAF.txt"),
        path("${patient}_${tumor_sample}_${normal_sample}_Germline_LogR.txt"),
        path("${patient}_${tumor_sample}_${normal_sample}_Tumor_BAF_rawBAF.txt"),
        path("${patient}_${tumor_sample}_${normal_sample}_Tumor_BAF.txt"),
        path("${patient}_${tumor_sample}_${normal_sample}_Tumor_LogR.txt")
    )

    script:
    """
    conda run -n CNAlignR Rscript /home/alg2264/repos/CNAlignR/scripts/wes_preprocess.R \
        --patient ${patient} \
        --tumor_name ${tumor_sample} \
        --tumor_bam ${tumor_bam} \
        --normal_name ${normal_sample} \
        --normal_bam ${normal_bam} \
        --sex ${sex} \
        --genome ${build} \
        --target_bed ${bed_file} \
        --allelecounter_exe ${allelecounter_exe} \
        --alleles_prefix ${alleles_prefix} \
        --loci_prefix ${loci_prefix} \
    """
}

/*
 * Get data object for CNalign
 */
process GET_CNALIGN_OBJ {
    tag "$patient"
    cpus 8
    memory '60GB'
    time '4h'
    executor 'slurm'
    queue 'short'
    publishDir params.ascat_dir, mode: 'copy'

    input:
    path all_allelecounter_files
    val normal_sample
    val patient
    val sex
    val build
    path GCcontentfile
    path replictimingfile

    output:
    //path "${patient}_CNalign_obj.rds"
    path "${patient}_CNalign_obj_mpcf.rds"
    //path "${patient}_CNalign_obj_mpcf_hisens.rds"

    script:
    """
    echo "Generating CNalign data object ..."
    conda run -n CNAlignR /home/alg2264/miniconda3/envs/CNAlignR/bin/Rscript /home/alg2264/repos/CNAlignR/scripts/wes_getinput.R \
        --patient ${patient} \
        --normal_sample ${normal_sample} \
        --sex ${sex} \
        --build ${build} \
        --GCcontentfile ${GCcontentfile} \
        --replictimingfile ${replictimingfile} \
        --obj_file "${patient}_CNalign_obj.rds"
    """
}

/*
 * Run snp-pileup for all samples
 */
process SNP_PILEUP {
    tag "$patient"
    cpus 1
    memory '32GB'
    time '8h'
    executor 'slurm'
    queue 'short'
    publishDir params.ascat_dir, mode: 'copy'

    input:
    path polymorphic_sites
    path polymorphic_sites_tbi
    path all_bams
    path all_bam_indices
    val patient

    output:
    path "${patient}.pileup.gz"

    script:
    def bam_line = all_bams.collect { bam -> "${bam}" }.join(' ')
    def n_bams   = all_bams.size()
    def r_arg    = (['25'] * n_bams).join(',')
    """
    conda run -n snp-pileup snp-pileup ${polymorphic_sites} ${patient}.pileup.gz \
        -q 10 -Q 20 -P 100 -d 4000 -g -r ${r_arg} ${bam_line}
    """
}

/*
 * Workflow
 */
workflow {

    // reference genome inputs
    ref_fasta = file(params.ref_fasta)
    ref_amb   = file(params.ref_amb)
    ref_ann   = file(params.ref_ann)
    ref_bwt   = file(params.ref_bwt)
    ref_fai   = file(params.ref_fai)
    ref_pac   = file(params.ref_pac)
    ref_sa    = file(params.ref_sa)
    ref_dict  = file(params.ref_dict)
    ref_files = tuple(ref_fasta, ref_amb, ref_ann, ref_bwt, ref_fai, ref_pac, ref_sa, ref_dict)

    // polymorphic sites
    polymorphic_sites     = file(params.polymorphic_sites)
    polymorphic_sites_tbi = file(params.polymorphic_sites_tbi)
    polymorphic_sites_files = tuple(polymorphic_sites, polymorphic_sites_tbi)

    // germline resources
    germline_resource     = file(params.germline_resource)
    germline_resource_tbi = file(params.germline_resource_tbi)
    germline_resource_files = tuple(germline_resource, germline_resource_tbi)

    // panel of normals
    panel_of_normals     = file(params.panel_of_normals)
    panel_of_normals_idx = file(params.panel_of_normals_idx)
    panel_of_normals_files = tuple(panel_of_normals, panel_of_normals_idx)

    // channel of genomic chunks
    genome_chunk_ch = Channel.fromPath(params.genome_chunks)
        .splitCsv(header: false, sep: ",", strip: true)
        .map { row -> tuple(row[0], row[1], row[2], row[3]) }

    // =============================
    // BAM preprocessing
    // =============================

    // make all expected output directories
    make_dirs_ch = MAKE_DIRS(
        params.bam_dir, params.maf_dir, params.mosdepth_dir,
        params.mtbam_dir, params.haplocheck_dir, params.ascat_dir,
        params.tmp_dir, params.mutect_dir, params.fastqc_dir
    )

    // Convert input BAMs to paired-end FASTQ
    bam_to_fastq_output = BAM_TO_FASTQ(sample_ch, make_dirs_ch)

    // FastQC (pre-trimming QC) + adapter trimming in one process.
    // Raw FASTQs are deleted inside the script once both steps are complete.
    fastqc_and_trim_output = FASTQC_AND_TRIM(bam_to_fastq_output)

    // BWA-MEM alignment (deletes trimmed FASTQs inside script)
    bwa_mem_output = BWA_MEM(fastqc_and_trim_output.trimmed, ref_files)

    // Sort SAM -> BAM (deletes raw SAM inside script)
    sortsam_output = SAMTOOLS_SORT(bwa_mem_output, ref_files)

    // MarkDuplicates (deletes sorted BAM/BAI inside script)
    markdup_output = GATK_MARKDUP(sortsam_output, ref_files)

    // Drop the metrics file from the tuple — downstream steps expect (sample, bam, bai)
    preprocessed_ch = markdup_output.map { sample, bam, bai, metrics -> tuple(sample, bam, bai) }

    // Extract BAM/index channels for multi-sample steps
    preprocessed_bam_ch       = preprocessed_ch.map { it[1] }
    preprocessed_bam_index_ch = preprocessed_ch.map { it[2] }

    // Collect all BAMs/indices for multi-sample steps
    all_bams_ch        = preprocessed_bam_ch.collect()
    all_bam_indices_ch = preprocessed_bam_index_ch.collect()

    // =============================
    // Variant calling & filtering
    // =============================

    mutect_output = GATK_MUTECT2(
        genome_chunk_ch, params.targets_bed,
        all_bams_ch, all_bam_indices_ch,
        params.patient, params.normal_sample, params.mutect_dir,
        polymorphic_sites_files, germline_resource_files, panel_of_normals_files,
        ref_files
    )

    region_vcf_ch   = mutect_output.map { it[2] }
    region_vcf_tbi_ch = mutect_output.map { it[3] }
    region_stats_ch = mutect_output.map { it[4] }
    region_f1r2_ch  = mutect_output.map { it[5] }

    all_vcf_ch      = region_vcf_ch.collect()
    all_vcf_tbi_ch  = region_vcf_tbi_ch.collect()
    all_stats_ch    = region_stats_ch.collect()
    all_f1r2_ch     = region_f1r2_ch.collect()

    // MERGE_REGIONS deletes per-region files inside script
    mergeregions_output = MERGE_REGIONS(
        params.patient, all_vcf_ch, all_vcf_tbi_ch, all_stats_ch, all_f1r2_ch
    )

    filter_calls_output = FILTER_MUTECT_CALLS(params.patient, mergeregions_output, ref_files)
    filtered_vcf_ch     = filter_calls_output.map { it[2] }
    filtered_vcf_tbi_ch = filter_calls_output.map { it[3] }

    // VCF2MAF (uses original sample_ch which now carries bam_dir instead of fq_dir)
    vcf2maf_output = VCF2MAF(sample_ch, filtered_vcf_ch, filtered_vcf_tbi_ch)

    // Mosdepth (QC)
    mosdepth_output = MOSDEPTH(preprocessed_ch, params.targets_bed)

    // Slice mtDNA
    slice_mtdna_output = SLICE_MTDNA(preprocessed_ch, params.mt_label)
    mtbam_ch       = slice_mtdna_output.map { it[1] }
    mtbam_index_ch = slice_mtdna_output.map { it[2] }
    all_mtbam_ch       = mtbam_ch.collect()
    all_mtbam_index_ch = mtbam_index_ch.collect()

    // CNA data prep (tumor/normal pairs)
    tumor_input_ch  = preprocessed_ch.filter { sample, bam, bai -> sample != params.normal_sample }
    normal_input_ch = preprocessed_ch.filter { sample, bam, bai -> sample == params.normal_sample }
    prep_input_ch   = tumor_input_ch.combine(normal_input_ch)

    prep_data_output = PREP_CNA_DATA(
        prep_input_ch, params.patient, params.sex, params.build,
        params.targets_bed, params.allelecounter_exe,
        params.alleles_prefix, params.loci_prefix
    )

    all_allelecounter_files_ch = prep_data_output.collect()

    cnalign_output = GET_CNALIGN_OBJ(
        all_allelecounter_files_ch, params.normal_sample,
        params.patient, params.sex, params.build,
        params.gccontentfile, params.replictimingfile
    )

    // SNP pileup
    //snp_pileup_output = SNP_PILEUP(
    //    params.polymorphic_sites, params.polymorphic_sites_tbi,
    //    all_bams_ch, all_bam_indices_ch, params.patient
    //)
}


