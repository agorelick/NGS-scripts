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
def sample_ch = Channel.from( rows.collect { [it.sample, it.input_mbam, it.input_mbam_index, it.order] } )

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
  tuple val(sample), path(input_mbam), path(input_mbam_index), val(sample_order)
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
    // single joined tuple: [sample, sample_order, raw_sam, reverted_bam, reverted_bam_sbi]
    tuple val(sample), val(sample_order), path(sample_raw_sam), path(reverted_mbam), path(reverted_mbam_sbi)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), val(sample_order), path("${sample}_merged.bam"), path("${sample}_merged.bam.bai")

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
    // single joined tuple: [sample, sample_order, shifted_raw_sam, reverted_bam, reverted_bam_sbi]
    tuple val(sample), val(sample_order), path(sample_raw_sam_shifted), path(reverted_mbam), path(reverted_mbam_sbi)
    tuple path(ref_shifted_fasta), path(ref_shifted_amb), path(ref_shifted_ann), path(ref_shifted_bwt), path(ref_shifted_fai), path(ref_shifted_pac), path(ref_shifted_sa), path(ref_shifted_dict)

    output:
    tuple val(sample), val(sample_order), path("${sample}_shifted_merged.bam"), path("${sample}_shifted_merged.bam.bai")

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

 
// 
// /*
//  * Run GATK Mutect2 for multi-sample tumor/normal variant calling
//  */
// process GATK_MUTECT2 {
// 
//     tag "$region"
//     cpus 18
//     memory '32GB'
//     time '24h'  // 2h is usually sufficient
//     executor 'slurm'
//     queue 'medium'
// 
//     publishDir params.mutect_dir, mode: 'copy'
// 
//     input:
//     tuple val(chr), val(start), val(end), val(region)               // split genome regions into equal sized chunks for parallelization
//     path bed_file
//     path all_bams
//     path all_bam_indices
//     val patient                                                     // patient ID
//     val normal_sample                                               // SM value for the normal sample
//     path output_dir                                                 // location for output files
//     tuple path(polymorphic_sites), path(polymorphic_sites_tbi)      // polymorphic sites
//     tuple path(germline_resource), path(germline_resource_tbi)      // germline resource
//     tuple path(panel_of_normals), path(panel_of_normals_idx)        // panel of normals
//     tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)  // ref genome files
// 
//     output:
//     tuple val(region), path("regions_${region}.bed"), path("${patient}_${region}.vcf.gz"), path("${patient}_${region}.vcf.gz.tbi"), path("${patient}_${region}.vcf.gz.stats"), path("${patient}_${region}.f1r2.tar.gz")
// 
//     script:
//     def bams_line = all_bams.collect { bam -> "-I ${bam}" }.join(' ')
//     """
// 
//     # subset the bed file for regions within the specified range
//     module load bedtools/2.31.0
//     echo -e "${chr}\t${start}\t${end}" | bedtools intersect -a ${bed_file} -b - > regions_${region}.bed
// 
//     # run mutect2 for this region
//     conda run -n gatk_4.6.1.0 gatk Mutect2 -R $ref_fasta \
//         $bams_line \
//         -normal $normal_sample \
//         -L regions_${region}.bed \
//         --f1r2-tar-gz ${patient}_${region}.f1r2.tar.gz \
//         --native-pair-hmm-threads 16 \
//         -O ${patient}_${region}.vcf.gz \
//         --germline-resource $germline_resource \
//         --panel-of-normals $panel_of_normals \
//     """
// }
//  
// /*
//  * filter mutect calls
//  */
// process FILTER_MUTECT_CALLS { 
// 
//     tag "$patient"
//     cpus 8
//     memory '16GB'
//     time '30m'
//     executor 'slurm'
//     queue 'short'
// 
//     publishDir params.mutect_dir, mode: 'copy'
// 
//     input:
//     val patient
//     tuple path(raw_vcf), path(raw_vcf_tbi), path(raw_vcf_stats), path(raw_artifact_priors)
//     tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)  // ref genome files
// 
//     output:
//     tuple path("${patient}_unfiltered_norm.vcf.gz"), path("${patient}_unfiltered_norm.vcf.gz.tbi"), path("${patient}_filtered.vcf.gz"), path("${patient}_filtered.vcf.gz.tbi")
// 
//     script:
//     """
// 
//     module load bcftools/1.21
// 
//     # FilterMutectCalls
//     conda run -n gatk_4.6.1.0 gatk FilterMutectCalls -R $ref_fasta -V $raw_vcf --orientation-bias-artifact-priors $raw_artifact_priors -O ${patient}_unfiltered.vcf.gz
// 
//     # IndexFeatureFile
//     conda run -n gatk_4.6.1.0 gatk IndexFeatureFile -I ${patient}_unfiltered.vcf.gz 
// 
//     # normalize VCF to split multi-allelic sites
//     bcftools norm --multiallelics -both --fasta-ref $ref_fasta ${patient}_unfiltered.vcf.gz | bcftools view -I -O z -o ${patient}_unfiltered_norm.vcf.gz -
// 
//     # indexing normalized VCF
//     conda run -n gatk_4.6.1.0 gatk IndexFeatureFile -I ${patient}_unfiltered_norm.vcf.gz 
// 
//     # filtering mutations
//     bcftools view -i "FILTER='PASS'" ${patient}_unfiltered_norm.vcf.gz | bcftools view -I -O z -o ${patient}_filtered.vcf.gz -
// 
//     # indexing filtered VCF
//     conda run -n gatk_4.6.1.0 gatk IndexFeatureFile -I ${patient}_filtered.vcf.gz 
//     """
// }
// 
// /*
//  * Run VCF2MAF on filtered VCF file for each sample
//  */
// process VCF2MAF {
// 
//     tag "$sample"
//     cpus 1
//     memory '16GB'
//     time '1h'
//     executor 'slurm'
//     queue 'short'
// 
//     publishDir params.maf_dir, mode: 'copy'
// 
// 
//     input:
//     tuple val(sample), val(fq_prefix), val(fq_dir), val(sample_order)
//     path filtered_vcf
//     path filtered_vcf_tbi
// 
//     output:
//     path "${sample}.maf"
// 
//     script:
//     """
//     ## subset the multi-sample VCF for this sample
//     module load bcftools/1.21
//     bcftools view $filtered_vcf -s $sample > ${sample}.vcf
// 
//     ## run vcfmaf to get a sample-specific maf
//     conda run -n vep perl /home/alg2264/repos/vcf2maf/vcf2maf.pl --input-vcf ${sample}.vcf --output-maf ${sample}.maf --tumor-id ${sample} --remap-chain /home/alg2264/repos/vcf2maf/data/hg38_to_GRCh38.chain
//     """
// }
// 


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
//     // panel of normals
//     panel_of_normals = file(params.panel_of_normals)
//     panel_of_normals_idx = file(params.panel_of_normals_idx)
//     panel_of_normals_files = tuple(panel_of_normals, panel_of_normals_idx)
// 
//     // channel of genomic chunks
//     genome_chunk_ch = Channel.fromPath(params.genome_chunks)
//                     .splitCsv(header: false, sep: ",", strip: true)
//                     .map { row -> tuple(row[0], row[1], row[2], row[3]) }
// 
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

    // BWA-MEM alignment (shifted)
    bwa_mem_shifted_output = BWA_MEM_SHIFTED(samtofastq_output, ref_shifted_files)

    // merge the bwa-mem output with revert-sam output by the sample name
    merge_input = bwa_mem_output
        .join(revert_sam_output, by: [0, 1])

    // merge the bwa-mem output with revert-sam output by the sample name
    merge_shifted_input = bwa_mem_shifted_output
        .join(revert_sam_output, by: [0, 1])

    // MergeBamAlignment
    merge_bam_output = MERGE_BAM_ALIGNMENTS(merge_input, ref_files)

    // (shifted) MergeBamAlignment
    merge_bam_shifted_output = MERGE_BAM_ALIGNMENTS_SHIFTED(merge_shifted_input, ref_shifted_files)





    // Run SortSam on valid samples
    //sortsam_output = SAMTOOLS_SORT(bwa_mem_output, ref_files)

    // Run MarkDuplicates on sorted bam files
    //markdup_output = GATK_MARKDUP(sortsam_output, ref_files)

    // Extracting each element into separate channels
    //preprocessed_sample_ch = applybqsr_output.map { it[0] }
    //preprocessed_bam_ch = applybqsr_output.map { it[1] }
    //preprocessed_bam_index_ch = applybqsr_output.map { it[2] }

    // Collect bams/indices
    //all_bams_ch = preprocessed_bam_ch.collect()
    //all_bam_indices_ch = preprocessed_bam_index_ch.collect()

    // =============================
    // variant calling, filtering
    // =============================

    // call mutations in multi-sample paired T/N mode
    //mutect_output = GATK_MUTECT2(genome_chunk_ch, params.targets_bed, all_bams_ch, all_bam_indices_ch, params.patient, params.normal_sample, params.mutect_dir, polymorphic_sites_files, germline_resource_files, panel_of_normals_files, ref_files)

    // Extracting each element into separate channels
    //region_ch = mutect_output.map { it[0] }
    //region_bed_ch = mutect_output.map { it[1] }
    //region_vcf_ch = mutect_output.map { it[2] }
    //region_vcf_tbi_ch = mutect_output.map { it[3] }
    //region_stats_ch = mutect_output.map { it[4] }
    //region_f1r2_ch = mutect_output.map { it[5] }

    // collect the bam files so that we can do multi-sample variant calling
    //all_vcf_ch = region_vcf_ch.collect()
    //all_vcf_tbi_ch = region_vcf_tbi_ch.collect()
    //all_stats_ch = region_stats_ch.collect()
    //all_f1r2_ch = region_f1r2_ch.collect()

    // merge results from mutect for each genomic chunk into a single file
    //mergeregions_output = MERGE_REGIONS(params.patient, all_vcf_ch, all_vcf_tbi_ch, all_stats_ch, all_f1r2_ch)

    // filter mutect calls
    //filter_calls_output = FILTER_MUTECT_CALLS(params.patient, mergeregions_output, ref_files)
    //filtered_vcf_ch = filter_calls_output.map { it[2] }
    //filtered_vcf_tbi_ch = filter_calls_output.map { it[3] }

    // VCF2MAF
    //vcf2maf_output = VCF2MAF(sample_ch, filtered_vcf_ch, filtered_vcf_tbi_ch)

}



