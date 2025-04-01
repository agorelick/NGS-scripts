/*
 * Pipeline parameters
 */

// define variables for this specific patient
params.sample = ["5A", "5B", "5C", "5D"]
params.sample_index = [0, 1, 2, 3]
params.normal_sample = '5A'
params.patient = 'LM31'

// Other files and directories
params.fq_dir = "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1136465519_WES_LM31_34_42/00_fastq"
params.bam_dir = "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1136465519_WES_LM31_34_42/bams"
params.output_dir = '/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1136465519_WES_LM31_34_42/output'
params.maf_dir = '/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1136465519_WES_LM31_34_42/mafs'
params.tmp_dir = "/n/scratch/users/a/alg2264"

// Genome reference files
params.ref_fasta = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
params.ref_amb = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.amb"
params.ref_ann = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.ann"
params.ref_bwt = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.bwt"
params.ref_fai = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.fai"
params.ref_pac = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.pac"
params.ref_sa = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa.sa"
params.ref_dict = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.dict"

// additional reference files with index files
params.polymorphic_sites = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz"
params.polymorphic_sites_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz.tbi"
params.germline_resource = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz"
params.germline_resource_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz.tbi"
params.panel_of_normals = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf"
params.panel_of_normals_idx = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf.idx"
params.targets_bed = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed"
params.genome_chunks = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38_50Mbchunks.csv"

/*
 * Trim Illumina Universal Adapters
 */
process TRIM_ADAPTERS {

    tag "$sample"
    cpus 8
    memory '16GB'
    time '4h'
    executor 'slurm'
    queue 'short'

    input:
    val sample
    val sample_index
    val fq_dir

    output:
    tuple val(sample), val(sample_index), path("${sample}_trimmed_R1.fastq.gz"), path("${sample}_trimmed_R2.fastq.gz")

    script:
    """
    module load gcc/9.2.0 python/3.8.12 cutadapt/4.1

    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --cores=8 \\
        -o ${sample}_trimmed_R1.fastq.gz -p ${sample}_trimmed_R2.fastq.gz \\
        ${fq_dir}/${sample}_R1_001.fastq.gz ${fq_dir}/${sample}_R2_001.fastq.gz
    """
}

/*
 * Run BWA-MEM
 */
process BWA_MEM {

    tag "$sample"
    cpus 8
    memory '16GB'
    time '8h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), val(sample_index), path(trimmed_fq1), path(trimmed_fq2)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_raw.sam")

    script:
    """
    module load gcc/6.2.0 bwa/0.7.15

    bwa mem -M -t 8 -R '@RG\\tID:${sample_index}\\tSM:${sample}\\tPL:Illumina' ${ref_fasta} ${trimmed_fq1} ${trimmed_fq2} > ${sample}_raw.sam

    """
}

/*
 * Run GATK MarkDuplicates
 */
process GATK_MARKDUP {

    tag "$sample"
    cpus 20
    memory '16GB'
    time '12h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), path(raw_sam)
    path bam_dir
    path tmp_dir
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_marked_dup.bam"), path("${sample}_marked_dup.bam.bai")

    script:
    """
    module load gcc/6.2.0 gatk/4.1.9.0

    gatk MarkDuplicatesSpark --input ${raw_sam} --output ${sample}_marked_dup.bam --tmp-dir ${tmp_dir} --reference ${ref_fasta} -M ${bam_dir}/${sample}_marked_dup_metrics.txt --conf 'spark.executor.cores=20'
    """
}

/*
 * Run GATK BaseRecalibrator
 */
process GATK_BASERECAL {

    tag "$sample"
    cpus 16
    memory '16B'
    time '8h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), path(markdup_bam), path(markdup_bam_index)
    tuple path(polymorphic_sites), path(polymorphic_sites_tbi)    
    path tmp_dir
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path(markdup_bam), path(markdup_bam_index), path("${sample}_recal_data.table")

    script:
    """
    module load gcc/6.2.0 gatk/4.1.9.0

    gatk BaseRecalibrator -I ${markdup_bam} -R ${ref_fasta} --known-sites ${polymorphic_sites} -O ${sample}_recal_data.table --tmp-dir ${tmp_dir}
    """
}

/*
 * Run GATK ApplyBQSR
 */
process GATK_APPLYBQSR {

    tag "$sample"
    cpus 16
    memory '16GB'
    time '8h'
    executor 'slurm'
    queue 'short'

    publishDir params.bam_dir, mode: 'copy'

    input:
    tuple val(sample), path(markdup_bam), path(markdup_bam_index), path(recal_data_table)
    path bam_dir
    path tmp_dir
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bai")

    script:
    """
    module load gcc/6.2.0 gatk/4.1.9.0

    gatk ApplyBQSR -R ${ref_fasta} -I ${markdup_bam} --bqsr-recal-file ${recal_data_table} -O ${sample}.bam --tmp-dir ${tmp_dir}

    """
}


/*
 * Run GATK Mutect2 for multi-sample tumor/normal variant calling
 */
process GATK_MUTECT2 {

    tag "$region"
    cpus 16
    memory '16GB'
    time '2h'
    executor 'slurm'
    queue 'short'

    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(chr), val(start), val(end), val(region)               // split genome regions into equal sized chunks for parallelization
    path bed_file
    path all_bams
    path all_bam_indices
    val patient                                                     // patient ID
    val normal_sample                                               // SM value for the normal sample
    path output_dir                                                 // location for output files
    tuple path(polymorphic_sites), path(polymorphic_sites_tbi)      // polymorphic sites
    tuple path(germline_resource), path(germline_resource_tbi)      // germline resource
    tuple path(panel_of_normals), path(panel_of_normals_idx)        // panel of normals
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)  // ref genome files

    output:
    tuple val(region), path("regions_${region}.bed"), path("${patient}_${region}.vcf.gz"), path("${patient}_${region}.vcf.gz.tbi"), path("${patient}_${region}.vcf.gz.stats"), path("${patient}_${region}.f1r2.tar.gz")


    script:
    def bams_line = all_bams.collect { bam -> "-I ${bam}" }.join(' ')
    """

    # subset the bed file for regions within the specified range
    module load gcc/9.2.0 bedtools/2.30.0

    echo -e "${chr}\t${start}\t${end}" | bedtools intersect -a ${bed_file} -b - > regions_${region}.bed

    # run mutect2 for this region
    module load gcc/6.2.0 gatk/4.1.9.0
    gatk Mutect2 -R $ref_fasta \
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
 */
process MERGE_REGIONS { 

    tag "$patient"
    cpus 8
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    input:
    val patient
    path tmp_dir
    path all_vcf
    path all_vcf_tbi
    path all_stats
    path all_f1r2

    publishDir params.output_dir, mode: 'copy'

    output:
    tuple path("${patient}_raw.vcf.gz"), path("${patient}_raw.vcf.gz.tbi"), path("${patient}_raw.vcf.gz.stats"), path("${patient}_raw.artifact-prior.tar.gz")

    script:
    def vcf_line = all_vcf.collect { vcf -> "${vcf}" }.join(' ')
    def stats_line = all_stats.collect { statsfile -> "--stats ${statsfile}" }.join(' ')
    def f1r2_line = all_f1r2.collect { f1r2file -> "-I ${f1r2file}" }.join(' ')

    """
    module load gcc/6.2.0 gatk/4.1.9.0 bcftools/1.13
 
    # combine the VCF file for each genomic chunk into a single VCF
    bcftools concat ${vcf_line} -a > ${patient}_raw_unsorted.vcf
    bcftools sort ${patient}_raw_unsorted.vcf -O z -o ${patient}_raw.vcf.gz
    rm ${patient}_raw_unsorted.vcf

    # index the combined VCF
    gatk IndexFeatureFile -I ${patient}_raw.vcf.gz 

    # combine the mutect-stats file for each chunk
    gatk MergeMutectStats ${stats_line} --output ${patient}_raw.vcf.gz.stats

    # learn read orientation bias
    gatk LearnReadOrientationModel ${f1r2_line} --output ${patient}_raw.artifact-prior.tar.gz

    """
}



/*
 * filter mutect calls
 */
process FILTER_MUTECT_CALLS { 

    tag "$patient"
    cpus 8
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.output_dir, mode: 'copy'

    input:
    val patient
    path tmp_dir
    tuple path(raw_vcf), path(raw_vcf_tbi), path(raw_vcf_stats), path(raw_artifact_priors)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)  // ref genome files

    output:
    tuple path("${patient}_unfiltered_norm.vcf.gz"), path("${patient}_unfiltered_norm.vcf.gz.tbi"), path("${patient}_filtered.vcf.gz"), path("${patient}_filtered.vcf.gz.tbi")

    script:
    """
    module load gcc/6.2.0 gatk/4.1.9.0 bcftools/1.13
 
    # FilterMutectCalls
    gatk FilterMutectCalls -R $ref_fasta -V $raw_vcf --orientation-bias-artifact-priors $raw_artifact_priors -O ${patient}_unfiltered.vcf.gz

    # IndexFeatureFile
    gatk IndexFeatureFile -I ${patient}_unfiltered.vcf.gz --tmp-dir $tmp_dir

    # normalize VCF to split multi-allelic sites
    bcftools norm --multiallelics -both --fasta-ref $ref_fasta ${patient}_unfiltered.vcf.gz | bcftools view -I -O z -o ${patient}_unfiltered_norm.vcf.gz -

    # indexing normalized VCF
    gatk IndexFeatureFile -I ${patient}_unfiltered_norm.vcf.gz --tmp-dir $tmp_dir

    # filtering mutations
    bcftools view -i "%FILTER='PASS'" ${patient}_unfiltered_norm.vcf.gz | bcftools view -I -O z -o ${patient}_filtered.vcf.gz -

    # indexing filtered VCF
    gatk IndexFeatureFile -I ${patient}_filtered.vcf.gz --tmp-dir $tmp_dir
    """
}



/*
 * Run VCF2MAF on filtered VCF file for each sample
 */
process VCF2MAF {

    tag "$sample"
    cpus 1
    memory '4GB'
    time '20m'
    executor 'slurm'
    queue 'short'

    publishDir params.maf_dir, mode: 'copy'


    input:
    val sample
    path filtered_vcf
    path filtered_vcf_tbi

    output:
    path "${sample}.maf"

    script:
    """
    module load gcc/9.2.0 bcftools/1.14 samtools/1.15.1

    ## subset the multi-sample VCF for this sample
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

    // polymorphic sites
    polymorphic_sites = file(params.polymorphic_sites)
    polymorphic_sites_tbi = file(params.polymorphic_sites_tbi)
    polymorphic_sites_files = tuple(polymorphic_sites, polymorphic_sites_tbi)    

    // germline resources
    germline_resource = file(params.germline_resource)
    germline_resource_tbi = file(params.germline_resource_tbi)
    germline_resource_files = tuple(germline_resource, germline_resource_tbi)

    // panel of normals
    panel_of_normals = file(params.panel_of_normals)
    panel_of_normals_idx = file(params.panel_of_normals_idx)
    panel_of_normals_files = tuple(panel_of_normals, panel_of_normals_idx)

    // Sample and index channels
    sample_ch = Channel.of(params.sample).flatten()
    sample_index_ch = Channel.of(params.sample_index).flatten()

    // channel of genomic chunks
    genome_chunk_ch = Channel.fromPath(params.genome_chunks)
                    .splitCsv(header: false, sep: ",", strip: true)
                    .map { row -> tuple(row[0], row[1], row[2], row[3]) }


    // =============================
    // bam preprocessing
    // =============================

    // Adapter trimming
    trim_adapt_output = TRIM_ADAPTERS(sample_ch, sample_index_ch, params.fq_dir)

    // BWA-MEM alignment
    bwa_mem_output = BWA_MEM(trim_adapt_output, ref_files)

    // Run MarkDuplicates only on valid samples
    markdup_output = GATK_MARKDUP(bwa_mem_output, params.bam_dir, params.tmp_dir, ref_files)

    // BaseRecalibrator
    baserecal_output = GATK_BASERECAL(markdup_output, polymorphic_sites_files, params.tmp_dir, ref_files)
    
    // ApplyBQSR
    applybqsr_output = GATK_APPLYBQSR(baserecal_output, params.bam_dir, params.tmp_dir, ref_files)

    // Extracting each element into separate channels
    preprocessed_sample_ch = applybqsr_output.map { it[0] }
    preprocessed_bam_ch = applybqsr_output.map { it[1] }
    preprocessed_bam_index_ch = applybqsr_output.map { it[2] }

    // Collect bams/indices
    all_bams_ch = preprocessed_bam_ch.collect()
    all_bam_indices_ch = preprocessed_bam_index_ch.collect()


    // =============================
    // variant calling, filtering
    // =============================

    // call mutations in multi-sample paired T/N mode
    mutect_output = GATK_MUTECT2(genome_chunk_ch, params.targets_bed, all_bams_ch, all_bam_indices_ch, params.patient, params.normal_sample, params.output_dir, polymorphic_sites_files, germline_resource_files, panel_of_normals_files, ref_files)

    // Extracting each element into separate channels
    region_ch = mutect_output.map { it[0] }
    region_bed_ch = mutect_output.map { it[1] }
    region_vcf_ch = mutect_output.map { it[2] }
    region_vcf_tbi_ch = mutect_output.map { it[3] }
    region_stats_ch = mutect_output.map { it[4] }
    region_f1r2_ch = mutect_output.map { it[5] }

    // collect the bam files so that we can do multi-sample variant calling
    all_vcf_ch = region_vcf_ch.collect()
    all_vcf_tbi_ch = region_vcf_tbi_ch.collect()
    all_stats_ch = region_stats_ch.collect()
    all_f1r2_ch = region_f1r2_ch.collect()

    // merge results from mutect for each genomic chunk into a single file
    mergeregions_output = MERGE_REGIONS(params.patient, params.tmp_dir, all_vcf_ch, all_vcf_tbi_ch, all_stats_ch, all_f1r2_ch)

    // filter mutect calls
    filter_calls_output = FILTER_MUTECT_CALLS(params.patient, params.tmp_dir, mergeregions_output, ref_files)
    filtered_vcf_ch = filter_calls_output.map { it[2] }
    filtered_vcf_tbi_ch = filter_calls_output.map { it[3] }

    // VCF2MAF
    vcf2maf_output = VCF2MAF(sample_ch, filtered_vcf_ch, filtered_vcf_tbi_ch)


}



