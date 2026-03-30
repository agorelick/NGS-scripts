if (!params.csv) {
    error "Missing required parameter: --csv"
}


// Load the CSV into a list of maps (1 per row)
// params.csv is defined from the nextflow run command via --csv <CSV_FILE>
def df = file(params.csv)
           .text
           .split('\n')
           .collect { it.trim() }
           .findAll { it }  // remove empty lines

def headers = df[0].split(',')  // column names
def rows = df[1..-1].collect { line -> 
    def fields = line.split(',')
    [ : ].withDefault { null }.tap { map ->
        headers.eachWithIndex { h, i ->
            if (i < fields.size()) {
                map[h] = fields[i]
            } else {
                map[h] = null  // or some default value
            }
        }
    }
}

// Just to test if parsing works
println "Parsed headers: ${headers}"
println "First row: ${rows[0]}"


// Create tuple channel: (sample, bam, index, order)
def sample_ch = Channel.from( rows.collect { [it.sample, it.order, it.input_mbam, it.input_mbam_index] } )

// Set scalar params using first row with non-null value
def first_full_row = rows.find { it.patient }  // or any other required field
params.patient        = first_full_row.patient
params.normal_sample = first_full_row.normal_sample
params.sex            = first_full_row.sex
params.nextflow_dir    = first_full_row.nextflow_dir
params.output_dir   = first_full_row.output_dir
params.mosdepth_dir = first_full_row.mosdepth

// dynamically generated parameters
params.realigned_mbam_dir = "${params.output_dir}/${params.patient}/realigned_mbams"
params.mutect_dir = "${params.output_dir}/${params.patient}/mutect"
params.maf_dir = "${params.output_dir}/${params.patient}/mafs"
params.tmp_dir = "${params.nextflow_dir}/${params.patient}/tmp_files"

// mtDNA reference genome files
params.ref_fasta = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.fasta"
params.ref_amb = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.fasta.amb"
params.ref_ann = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.fasta.ann"
params.ref_bwt = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.fasta.bwt"
params.ref_fai = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.fasta.fai"
params.ref_pac = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.fasta.pac"
params.ref_sa = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.fasta.sa"
params.ref_dict = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.dict"

// (shifted) mtDNA reference genome files
params.ref_shifted_fasta = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
params.ref_shifted_amb = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb"
params.ref_shifted_ann = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann"
params.ref_shifted_bwt = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt"
params.ref_shifted_fai = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai"
params.ref_shifted_pac = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac"
params.ref_shifted_sa = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa"
params.ref_shifted_dict = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict"

params.control_region_shifted_intervals = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/control_region_shifted.chrM.interval_list"
params.non_control_region_intervals = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/non_control_region.chrM.interval_list"
params.shift_back_chain = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/mtdna/ShiftBack.chain"


// additional reference files with index files
// params.polymorphic_sites = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz"
// params.polymorphic_sites_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz.tbi"
// params.germline_resource = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz"
// params.germline_resource_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz.tbi"
// params.panel_of_normals = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf"
// params.panel_of_normals_idx = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf.idx"
 
process MAKE_DIRS {
    tag "mkdirs"
    executor 'local'

    input:
    val realigned_mbam_dir
    val mutect_dir
    val maf_dir
    val tmp_dir

    output:
    path "mkdir_done.txt"

    script:
    """
    mkdir -p ${realigned_mbam_dir}
    mkdir -p ${mutect_dir}
    mkdir -p ${maf_dir}
    mkdir -p ${tmp_dir}
    touch mkdir_done.txt
    """
}


process CLEAN_CHRM_BAM_TO_FASTQ {

  tag "$sample"
  cpus 2
  memory '8GB'
  time '1h'
  executor 'slurm'
  queue 'short'

  input:
  tuple val(sample), val(sample_order), path(input_mbam), path(input_mbam_index)
  val make_dirs_ch

  output:
  tuple val(sample), val(sample_order), path("${sample}.paired.R1.fastq.gz"), path("${sample}.paired.R2.fastq.gz")

  script:
  """
  # Name-collate so mates are adjacent
  samtools collate -@ ${task.cpus} -O ${input_mbam} ${sample}.namecollated | \
    samtools fastq \
      -@ ${task.cpus} \
      -1 ${sample}.paired.R1.fastq.gz \
      -2 ${sample}.paired.R2.fastq.gz \
      -s ${sample}.singletons.fastq.gz \
      -0 /dev/null \
      -n \
      -
  """
}

process SAMTOFASTQ {

  tag "$sample"
  cpus 1
  memory '8GB'
  time '20m'
  executor 'slurm'
  queue 'short'

  input:
  tuple val(sample), val(sample_order), path(reverted_mbam)

  output:
  tuple val(sample), val(sample_order), path("${sample}_R1.fastq"), path("${sample}_R2.fastq")

  script:
  """
  conda run -n gatk_4.6.1.0 picard SamToFastq \
    -I ${reverted_mbam} \
    -F ${sample}_R1.fastq \
    -F2 ${sample}_R2.fastq
  """
}



process FASTQ_TO_UBAM {

  tag "$sample"
  cpus 2
  memory '8GB'
  time '1h'
  executor 'slurm'
  queue 'short'

  input:
  tuple val(sample), val(sample_order), path(r1_fastq), path(r2_fastq)

  output:
  tuple val(sample), val(sample_order), path("${sample}_mt_clean.ubam")

  script:
  """
  conda run -n gatk_4.6.1.0 picard FastqToSam \
    FASTQ=${r1_fastq} \
    FASTQ2=${r2_fastq} \
    OUTPUT=${sample}_mt_clean.ubam \
    SAMPLE_NAME=${sample} \
    READ_GROUP_NAME=${sample} \
    LIBRARY_NAME=${sample} \
    PLATFORM=ILLUMINA \
    TMP_DIR=${params.tmp_dir}
  """
}


/*
 * Run BWA-MEM
 */
process BWA_MEM {

    tag "$sample"
    cpus 8
    memory '16GB'
    time '20m'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), val(sample_order), path(fq1), path(fq2)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), val(sample_order), path("${sample}_raw.sam")

    script:
    """
    module load gcc/14.2.0 bwa/0.7.18

    bwa mem -M -t 8 -R '@RG\\tID:${sample_order}\\tSM:${sample}\\tPL:Illumina' ${ref_fasta} ${fq1} ${fq2} > ${sample}_raw.sam
    """
}


/*
 * Run (shifted) BWA-MEM
 */
process BWA_MEM_SHIFTED {

    tag "$sample"
    cpus 8
    memory '16GB'
    time '20m'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), val(sample_order), path(fq1), path(fq2)
    tuple path(ref_shifted_fasta), path(ref_shifted_amb), path(ref_shifted_ann), path(ref_shifted_bwt), path(ref_shifted_fai), path(ref_shifted_pac), path(ref_shifted_sa), path(ref_shifted_dict)

    output:
    tuple val(sample), val(sample_order), path("${sample}_shifted_raw.sam")

    script:
    """
    module load gcc/14.2.0 bwa/0.7.18
    bwa mem -M -t 8 -R '@RG\\tID:${sample_order}\\tSM:${sample}\\tPL:Illumina' ${ref_shifted_fasta} ${fq1} ${fq2} > ${sample}_shifted_raw.sam
    """
}



/*
 * MERGE_BAM_ALIGNMENTS
 */
process MERGE_BAM_ALIGNMENTS {

    tag "$sample"
    cpus 4
    memory '8GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.realigned_mbam_dir, mode: 'copy'

    input:
    // single joined tuple: [sample, sample_order, raw_sam, reverted_bam]
    tuple val(sample), val(sample_order), path(sample_raw_sam), path(reverted_mbam)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), val(sample_order), path("${sample}_merged.bam"), path("${sample}_merged.bai")

    script:
    """
    conda run -n gatk_4.6.1.0 picard MergeBamAlignment \\
        -ALIGNED ${sample_raw_sam} \\
        -UNMAPPED ${reverted_mbam} \\
        -R ${ref_fasta} \\
        -O ${sample}_merged.bam \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    """
}


/*
 * MERGE_BAM_ALIGNMENTS_SHIFTED
 */
process MERGE_BAM_ALIGNMENTS_SHIFTED {

    tag "$sample"
    cpus 4
    memory '8GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.realigned_mbam_dir, mode: 'copy'

    input:
    // single joined tuple: [sample, sample_order, shifted_raw_sam, reverted_bam]
    tuple val(sample), val(sample_order), path(sample_raw_sam_shifted), path(reverted_mbam)
    tuple path(ref_shifted_fasta), path(ref_shifted_amb), path(ref_shifted_ann), path(ref_shifted_bwt), path(ref_shifted_fai), path(ref_shifted_pac), path(ref_shifted_sa), path(ref_shifted_dict)

    output:
    tuple val(sample), val(sample_order), path("${sample}_shifted_merged.bam"), path("${sample}_shifted_merged.bai")

    script:
    """
    conda run -n gatk_4.6.1.0 picard MergeBamAlignment \\
        -ALIGNED ${sample_raw_sam_shifted} \\
        -UNMAPPED ${reverted_mbam} \\
        -R ${ref_shifted_fasta} \\
        -O ${sample}_shifted_merged.bam \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    """
}



/*
 * Fix read group tags in merged BAM to ensure all reads have an RG tag
 * that matches the @RG header entry. This is necessary because MergeBamAlignment
 * can carry forward stale RG tags from the original whole-genome alignment
 * that have no corresponding @RG header entry in the chrM-realigned BAM,
 * causing Mutect2 to fail with "null sample name" errors.
 */

process FIX_READ_GROUPS {

    tag "$sample"
    cpus 2
    memory '8GB'
    time '30m'
    executor 'slurm'
    queue 'short'
    publishDir params.realigned_mbam_dir, mode: 'copy'

    input:
    tuple val(sample), val(sample_order), path(merged_bam), path(merged_bam_bai)

    output:
    tuple val(sample), val(sample_order), path("${sample}_merged_fixed.bam"), path("${sample}_merged_fixed.bai")

    script:
    """
    conda run -n gatk_4.6.1.0 picard AddOrReplaceReadGroups \\
        -I ${merged_bam} \\
        -O ${sample}_merged_fixed.bam \\
        --RGID ${sample} \\
        --RGSM ${sample} \\
        --RGLB ${sample} \\
        --RGPL ILLUMINA \\
        --RGPU 1 \\
        --QUIET true \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    """
}

/*
 * Same RG fix for the shifted-reference merged BAM.
 */

process FIX_READ_GROUPS_SHIFTED {

    tag "$sample"
    cpus 2
    memory '8GB'
    time '30m'
    executor 'slurm'
    queue 'short'
    publishDir params.realigned_mbam_dir, mode: 'copy'

    input:
    tuple val(sample), val(sample_order), path(merged_bam_shifted), path(merged_bam_shifted_bai)

    output:
    tuple val(sample), val(sample_order), path("${sample}_shifted_merged_fixed.bam"), path("${sample}_shifted_merged_fixed.bai")

    script:
    """
    conda run -n gatk_4.6.1.0 picard AddOrReplaceReadGroups \\
        -I ${merged_bam_shifted} \\
        -O ${sample}_shifted_merged_fixed.bam \\
        --RGID ${sample} \\
        --RGSM ${sample} \\
        --RGLB ${sample} \\
        --RGPL ILLUMINA \\
        --RGPU 1 \\
        --QUIET true \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    """
}




/*
 * Run Mutect2 on the non-shifted (non-control region) BAM
 */

process MUTECT2 {

publishDir params.mutect_dir, mode: 'copy'

    tag "$patient"
    cpus 4
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    input:
    path all_bams
    path all_bam_indices
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)
    path non_control_region_intervals
    val patient

    output:
    tuple path("${patient}_nonshifted.vcf.gz"), path("${patient}_nonshifted.vcf.gz.tbi"), path("${patient}_nonshifted.vcf.gz.stats"), path("${patient}_nonshifted.f1r2.tar.gz")

    script:
    def bams_line = all_bams.collect { bam -> "-I ${bam}" }.join(' ')
    """
    conda run -n gatk_4.6.1.0 gatk Mutect2 \\
        -R ${ref_fasta} \\
        ${bams_line} \\
        --mitochondria-mode \\
        -L ${non_control_region_intervals} \\
        --annotation StrandBiasBySample \\
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \\
        --read-filter MateUnmappedAndUnmappedReadFilter \\
        --max-reads-per-alignment-start 75 \\
        --max-mnp-distance 0 \\
        --f1r2-tar-gz ${patient}_nonshifted.f1r2.tar.gz \\
        -O ${patient}_nonshifted.vcf.gz
    """
}



/*
 * Run Mutect2 on the shifted (control region) BAM, then liftover
 * the resulting VCF back to standard chrM coordinates.
 *
 * Mutect2 is run in --mitochondria-mode on the shifted-reference BAM,
 * restricted to the control region interval list (in shifted coordinates).
 * LiftoverVcf then maps variant positions back to standard chrM coords
 * using the shift-back chain file.
 */

process MUTECT2_SHIFTED {

    publishDir params.mutect_dir, mode: 'copy'

    tag "$patient"
    cpus 4
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    input:
    path all_shifted_bams
    path all_shifted_bam_indices
    tuple path(ref_shifted_fasta), path(ref_shifted_amb), path(ref_shifted_ann), path(ref_shifted_bwt), path(ref_shifted_fai), path(ref_shifted_pac), path(ref_shifted_sa), path(ref_shifted_dict)
    path control_region_shifted_intervals
    val patient

    output:
    tuple path("${patient}_shifted.vcf.gz"), path("${patient}_shifted.vcf.gz.tbi"), path("${patient}_shifted.vcf.gz.stats"), path("${patient}_shifted.f1r2.tar.gz")

    script:
    def shifted_bams_line = all_shifted_bams.collect { bam -> "-I ${bam}" }.join(' ')
    """
    conda run -n gatk_4.6.1.0 gatk Mutect2 \\
        -R ${ref_shifted_fasta} \\
        ${shifted_bams_line} \\
        --mitochondria-mode \\
        -L ${control_region_shifted_intervals} \\
        --annotation StrandBiasBySample \\
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \\
        --read-filter MateUnmappedAndUnmappedReadFilter \\
        --max-reads-per-alignment-start 75 \\
        --max-mnp-distance 0 \\
        --f1r2-tar-gz ${patient}_shifted.f1r2.tar.gz \\
        -O ${patient}_shifted.vcf.gz
    """
}

/*
 * Liftover the shifted Mutect2 VCF back to standard chrM coordinates.
 *
 * Uses the shift-back chain file generated by GATK ShiftFasta to map
 * variant positions from shifted chrM coordinates back to standard chrM.
 * Variants that cannot be lifted over are written to a separate reject file
 * for inspection.
 *
 * --RECOVER_SWAPPED_REF_ALT: rescues variants where ref/alt are swapped
 * after liftover due to strand orientation differences between the two
 * reference versions.
 */

process LIFTOVER_VCF {

    publishDir params.mutect_dir, mode: 'copy'

    tag "$patient"
    cpus 2
    memory '8GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    input:
    tuple path(shifted_vcf), path(shifted_vcf_tbi), path(shifted_vcf_stats), path(shifted_f1r2)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)
    path shift_back_chain
    val patient

    output:
    tuple path("${patient}_shifted_back.vcf.gz"), path("${patient}_shifted_back.vcf.gz.tbi"), path(shifted_vcf_stats), path(shifted_f1r2)

    script:
    """
    conda run -n gatk_4.6.1.0 picard LiftoverVcf \\
        -I ${shifted_vcf} \\
        -O ${patient}_shifted_back.vcf.gz \\
        --CHAIN ${shift_back_chain} \\
        --REJECT ${patient}_liftover_rejected.vcf.gz \\
        -R ${ref_fasta} \\
        --RECOVER_SWAPPED_REF_ALT true \\
        --VALIDATION_STRINGENCY LENIENT \\
        --CREATE_INDEX true
    """
}


/*
 * Merge the standard (non-control region) and shifted-back (control region)
 * VCFs into a single whole-chrM VCF, then apply FilterMutectCalls.
 *
 * MergeVcfs combines the two complementary callsets into one VCF spanning
 * the entire mitochondrial genome. The input VCFs must share the same
 * coordinate system (both in standard chrM coords after liftover).
 *
 * FilterMutectCalls applies mitochondria-specific filters:
 *   --mitochondria-mode      : sets appropriate filter thresholds for MT
 *   --stats                  : merged stats from both Mutect2 runs, used to
 *                              train the somatic/artifact model
 *   --ob-priors              : orientation bias artifact priors from
 *                              LearnReadOrientationModel (FFPE/OxoG)
 *   --min-allele-fraction    : hard VAF floor (0.01 = 1%)
 *   --max-alt-allele-count   : filter noisy multiallelic pileups
 *
 * Note: MergeStats and LearnReadOrientationModel must be run before this
 * process. Their outputs (merged stats, artifact priors) are passed in
 * alongside the two VCFs.
 */

process MERGE_VCFS_AND_FILTER {

    tag "$patient"
    cpus 2
    memory '8GB'
    time '1h'
    executor 'slurm'
    queue 'short'

    publishDir params.mutect_dir, mode: 'copy'

    input:
    tuple path(vcf), path(vcf_tbi), path(vcf_stats), path(f1r2)
    tuple path(vcf_shifted_back), path(vcf_shifted_back_tbi), path(vcf_shifted_back_stats), path(f1r2_shifted)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)
    val patient

    output:
    tuple path("${patient}_filtered.vcf.gz"), path("${patient}_filtered.vcf.gz.tbi")

    script:
    """
    # Step 1: Merge stats files from both Mutect2 runs
    conda run -n gatk_4.6.1.0 gatk MergeMutectStats \\
        --stats ${vcf_stats} \\
        --stats ${vcf_shifted_back_stats} \\
        -O ${patient}_merged.stats

    # Step 2: Learn orientation bias artifact priors from both f1r2 tarballs
    conda run -n gatk_4.6.1.0 gatk LearnReadOrientationModel \\
        -I ${f1r2} \\
        -I ${f1r2_shifted} \\
        -O ${patient}_artifact_priors.tar.gz

    # Step 3: Merge the two VCFs into a single whole-chrM VCF
    conda run -n gatk_4.6.1.0 picard MergeVcfs \\
        -I ${vcf} \\
        -I ${vcf_shifted_back} \\
        -O ${patient}_merged.vcf.gz \\
        --SEQUENCE_DICTIONARY ${ref_dict} \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT

    # Step 4: Filter the merged VCF
    conda run -n gatk_4.6.1.0 gatk FilterMutectCalls \\
        -V ${patient}_merged.vcf.gz \\
        -R ${ref_fasta} \\
        -O ${patient}_filtered.vcf.gz \\
        --stats ${patient}_merged.stats \\
        --mitochondria-mode \\
        --ob-priors ${patient}_artifact_priors.tar.gz \\
        --min-allele-fraction 0.01 \\
        --max-alt-allele-count 4 \\
        --create-output-variant-index true
    """
}




/*
 * Run VCF2MAF on filtered VCF file for each sample
 */
process VCF2MAF {

    tag "$sample"
    cpus 1
    memory '16GB'
    time '20m'
    executor 'slurm'
    queue 'short'

    publishDir params.maf_dir, mode: 'copy'

    input:
    tuple val(sample), val(fq_prefix), val(fq_dir), val(sample_order)
    tuple path(filtered_vcf), path(filtered_vcf_tbi)

    output:
    path "${sample}.maf"

    script:
    """
    ## subset the multi-sample VCF for this sample
    module load bcftools/1.21
    bcftools view $filtered_vcf -s $sample > ${sample}.vcf

    ## run vcfmaf to get a sample-specific maf
    conda run -n vep perl /home/alg2264/repos/vcf2maf/vcf2maf.pl --input-vcf ${sample}.vcf --output-maf ${sample}.maf --tumor-id ${sample} --remap-chain /home/alg2264/repos/vcf2maf/data/hg38_to_GRCh38.chain
    """
}



/*
 * Workflow
 */
workflow {

     // reference genome inputs
     ref_fasta = file(params.ref_fasta)
     ref_amb = file(params.ref_amb)
     ref_ann = file(params.ref_ann)
     ref_bwt = file(params.ref_bwt)
     ref_fai = file(params.ref_fai)
     ref_pac = file(params.ref_pac)
     ref_sa = file(params.ref_sa)
     ref_dict = file(params.ref_dict)
     ref_files = tuple(ref_fasta, ref_amb, ref_ann, ref_bwt, ref_fai, ref_pac, ref_sa, ref_dict)
 
     // reference genome inputs
     ref_shifted_fasta = file(params.ref_shifted_fasta)
     ref_shifted_amb = file(params.ref_shifted_amb)
     ref_shifted_ann = file(params.ref_shifted_ann)
     ref_shifted_bwt = file(params.ref_shifted_bwt)
     ref_shifted_fai = file(params.ref_shifted_fai)
     ref_shifted_pac = file(params.ref_shifted_pac)
     ref_shifted_sa = file(params.ref_shifted_sa)
     ref_shifted_dict = file(params.ref_shifted_dict)
     ref_shifted_files = tuple(ref_shifted_fasta, ref_shifted_amb, ref_shifted_ann, ref_shifted_bwt, ref_shifted_fai, ref_shifted_pac, ref_shifted_sa, ref_shifted_dict)

     control_region_shifted_intervals = file(params.control_region_shifted_intervals)
     non_control_region_intervals = file(params.non_control_region_intervals)

 
//     // polymorphic sites
//     polymorphic_sites = file(params.polymorphic_sites)
//     polymorphic_sites_tbi = file(params.polymorphic_sites_tbi)
//     polymorphic_sites_files = tuple(polymorphic_sites, polymorphic_sites_tbi)    
// 
//     // germline resources
//     germline_resource = file(params.germline_resource)
//     germline_resource_tbi = file(params.germline_resource_tbi)
//     germline_resource_files = tuple(germline_resource, germline_resource_tbi)
// 

    // =============================
    // bam preprocessing
    // =============================

    // run process to make all the expected directories for this patient
    make_dirs_ch = MAKE_DIRS(params.realigned_mbam_dir, params.mutect_dir, params.maf_dir, params.tmp_dir)

    clean_fastq_output = CLEAN_CHRM_BAM_TO_FASTQ(sample_ch, make_dirs_ch)
    ubam_output        = FASTQ_TO_UBAM(clean_fastq_output)
    samtofastq_output  = SAMTOFASTQ(ubam_output)

    // BWA-MEM alignment
    bwa_mem_output = BWA_MEM(samtofastq_output, ref_files)
    bwa_mem_shifted_output = BWA_MEM_SHIFTED(samtofastq_output, ref_shifted_files)

    // merge the bwa-mem output with revert-sam output by the sample name
    merge_input = bwa_mem_output
        .join(ubam_output, by: [0, 1])
    merge_shifted_input = bwa_mem_shifted_output
        .join(ubam_output, by: [0, 1])

    // MergeBamAlignment
    merge_bam_output = MERGE_BAM_ALIGNMENTS(merge_input, ref_files)
    merge_bam_shifted_output = MERGE_BAM_ALIGNMENTS_SHIFTED(merge_shifted_input, ref_shifted_files)

    // Fix stale RG tags before variant calling
    fix_rg_output = FIX_READ_GROUPS(merge_bam_output)
    fix_rg_shifted_output = FIX_READ_GROUPS_SHIFTED(merge_bam_shifted_output)


    // Collect bams/indices
    sample_ch = fix_rg_output.map { it[0] }
    sample_order_ch = fix_rg_output.map { it[1] }
    bam_ch = fix_rg_output.map { it[2] }
    bam_index_ch = fix_rg_output.map { it[3] }
    all_bams_ch = bam_ch.collect()
    all_bam_indices_ch = bam_index_ch.collect()
    
    
    // Collect shifted bams/indices
    shifted_sample_ch = fix_rg_shifted_output.map { it[0] }
    shifted_sample_order_ch = fix_rg_shifted_output.map { it[1] }
    shifted_bam_ch = fix_rg_shifted_output.map { it[2] }
    shifted_bam_index_ch = fix_rg_shifted_output.map { it[3] }
    all_shifted_bams_ch = shifted_bam_ch.collect()
    all_shifted_bam_indices_ch = shifted_bam_index_ch.collect()


    // call variants in the non-shifted VCF
    mutect2_output = MUTECT2(all_bams_ch, all_bam_indices_ch, ref_files, non_control_region_intervals, params.patient)
    mutect2_shifted_output = MUTECT2_SHIFTED(all_shifted_bams_ch, all_shifted_bam_indices_ch, ref_shifted_files, control_region_shifted_intervals, params.patient)
    
    // lift over the shifted VCF 
    liftover_output = LIFTOVER_VCF(mutect2_shifted_output, ref_files, params.shift_back_chain, params.patient)

    // merge standard and shifted VCFs, apply filtering
    filtered_output = MERGE_VCFS_AND_FILTER(mutect2_output, liftover_output, ref_files, params.patient)

    // VCF2MAF
    vcf2maf_output = VCF2MAF(sample_ch, filtered_output)

}



