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

def headers = df[0].split(',').collect { it.trim() }

def rows = df[1..-1].collect { line ->
    def fields = line.split(',').collect { it.trim() }
    def row = [:].withDefault { null }
    headers.eachWithIndex { h, i ->
        if (i < fields.size()) {
            row[h] = fields[i]
        }
    }
    return row
}

// Set scalar params using first row with non-null value
def first_full_row = rows.find { it.patient }  // or any other required field
params.nextflow_dir    = first_full_row.nextflow_dir
params.output_dir   = first_full_row.output_dir
params.patient   = first_full_row.patient

// dynamically generated parameters
params.maf_dir = "${params.output_dir}/${params.patient}/mafs"
params.mito_output_dir = "${params.output_dir}/${params.patient}/mtdna"
params.bam_dir = "${params.output_dir}/${params.patient}/bams"
params.mtbam_dir = "${params.output_dir}/${params.patient}/mtbams"
params.tmp_dir = "${params.nextflow_dir}/${params.patient}/tmp_files"

// Just to test if parsing works
println "[INFO] Loaded Nextflow params:"
params.each { k, v ->
    println "  ${k} = ${v}"
}

def sample_ch = Channel.from( rows.collect { [it.sample, it.bam, it.bam_index, it.barcode] } )

/*
 * Additional pipeline parameters (use for all WES data)
 */

// Genome reference files
params.build = "hg38"
params.mt_label = "chrM"
params.ref_fasta = "/home/alg2264/nextflow/snrna_multisample/ff01/genome.fa"
params.ref_fai = "/home/alg2264/nextflow/snrna_multisample/ff01/genome.fa.fai"
params.ref_dict = "/home/alg2264/nextflow/snrna_multisample/ff01/genome.dict"

// temporary, replace this once we add split bam process
params.split_bam_dir = "/n/data1/hms/genetics/naxerova/lab/alex/peritoneal_revision_WES/snRNAseq/NaxID_ff_01/split_bams"

// additional reference files with index files
params.polymorphic_sites = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz"
params.polymorphic_sites_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz.tbi"
params.germline_resource = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz"
params.germline_resource_tbi = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/gnomad.raw.sites.hg38/af-only-gnomad.hg38.vcf.gz.tbi"
params.panel_of_normals = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf"
params.panel_of_normals_idx = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/PoN/1000g_pon.hg38.vcf.idx"
params.targets_bed = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38.bed"
params.genome_chunks = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38_50Mbchunks.csv"
params.rediportal = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/REDIportal.bed.gz"

process MAKE_DIRS {
    tag "mkdirs"
    executor 'local'

    input:
    val output_dir
    val maf_dir
    val mtbam_dir
    val tmp_dir
    val bam_dir


    output:
    path "mkdir_done.txt"

    """
    mkdir -p ${output_dir}
    mkdir -p ${maf_dir}
    mkdir -p ${mtbam_dir}
    mkdir -p ${tmp_dir}
    mkdir -p ${bam_dir}

    touch mkdir_done.txt
    """
}


/*
 * Run GATK Split'N'Trim & Reassign Mapping Qualities
 */
process GATK_SPLITNTRIM {

    tag "$sample"
    cpus 8
    memory '20GB'
    time '4h' 
    executor 'slurm'
    queue 'short'

    publishDir params.bam_dir, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bam_index), path(barcode)
    val patient
    path bam_dir
    path split_bam_dir
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)  // ref genome files

    output:
    tuple val(sample), path("${sample}_tumor.split.bam"), path("${sample}_tumor.split.bai"), path("${sample}_normal.split.bam"), path("${sample}_normal.split.bai")

    script:
    """
    # subset the bed file for regions within the specified range
    module load gcc/6.2.0 gatk/4.1.9.0
 
    gatk SplitNCigarReads -R $ref_fasta -I ${split_bam_dir}/${sample}_tumor.bam -O ${sample}_tumor.split.bam
    gatk SplitNCigarReads -R $ref_fasta -I ${split_bam_dir}/${sample}_normal.bam -O ${sample}_normal.split.bam
    """
}


/*
 * Run GATK merge to create a single merged normal bam
 */
process GATK_MERGE_NORMAL {

    tag "$patient"
    cpus 8
    memory '20GB'
    time '2h'
    executor 'slurm'
    queue 'short'

    publishDir params.bam_dir, mode: 'copy'

    input:
    path all_normal_bams
    path all_normal_bam_indices
    val patient

    output:
    tuple path("${patient}_normal_merged.bam"), path("${patient}_normal_merged.bam.bai")

    script:
    def bams_line = all_normal_bams.collect { bam -> "${bam}" }.join(' ')
    """

    # subset the bed file for regions within the specified range
    module load gcc/9.2.0 samtools/1.15.1

    samtools merge -r -o ${patient}_normal_merged_tmp.bam ${bams_line}
    samtools view -H ${patient}_normal_merged_tmp.bam | sed "s/SM:[^\t]*/SM:Nall/g" | samtools reheader - ${patient}_normal_merged_tmp.bam > ${patient}_normal_merged.bam
    samtools index ${patient}_normal_merged.bam
    """
}




/*
 * Run GATK Mutect2 for multi-sample tumor/normal variant calling
 */
process GATK_MUTECT2 {

    tag "$region"
    cpus 18
    memory '32GB'
    time '2h'  // 2h is usually sufficient
    executor 'slurm'
    queue 'short'

    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(chr), val(start), val(end), val(region)               // split genome regions into equal sized chunks for parallelization
    path bed_file
    path all_bams
    path all_bam_indices
    tuple path(normal_bam), path(normal_bam_index)
    val patient                                                     // patient ID
    path output_dir                                                 // location for output files
    tuple path(polymorphic_sites), path(polymorphic_sites_tbi)      // polymorphic sites
    tuple path(germline_resource), path(germline_resource_tbi)      // germline resource
    tuple path(panel_of_normals), path(panel_of_normals_idx)        // panel of normals
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)  // ref genome files

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
        $bams_line -I ${normal_bam} \
        -normal Nall \
        -L regions_${region}.bed \
        --f1r2-tar-gz ${patient}_${region}.f1r2.tar.gz \
        --native-pair-hmm-threads 16 \
        -O ${patient}_${region}.vcf.gz \
        --downsampling-stride 20 \
        --max-reads-per-alignment-start 6 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --germline-resource $germline_resource \
        --panel-of-normals $panel_of_normals \
    """
}


/*
 * Run GATK Mutect2 (mitochondrial mode) for each sample
 */
process GATK_MUTECT2_MITO {

    tag "$sample"
    cpus 18
    memory '32GB'
    time '45m'  // 2h is usually sufficient
    executor 'slurm'
    queue 'short'

    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(sample), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index)
    val patient
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)            // ref genome files
    path output_dir                                                 // location for output files
    tuple path(polymorphic_sites), path(polymorphic_sites_tbi)      // polymorphic sites
    tuple path(germline_resource), path(germline_resource_tbi)      // germline resource
    tuple path(panel_of_normals), path(panel_of_normals_idx)        // panel of normals

    output:
    tuple val(sample), path("${patient}_${sample}_chrM.vcf.gz"), path("${patient}_${sample}_chrM.vcf.gz.tbi"), path("${patient}_${sample}_chrM.vcf.gz.stats"), path("${patient}_${sample}_chrM.f1r2.tar.gz")

    script:
    """

    # run mutect2 for this region
    module load gcc/6.2.0 gatk/4.1.9.0

    gatk Mutect2 -R $ref_fasta \
        -L chrM \
        --mitochondria-mode true \
        -I ${tumor_bam} \
        -O ${patient}_${sample}_chrM.vcf.gz \
        --disable-read-filter MappingQualityAvailableReadFilter \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --f1r2-tar-gz ${patient}_${sample}_chrM.f1r2.tar.gz \
        --native-pair-hmm-threads 16
    """
}


/*
 * Run GATK Mutect2 (mitochondrial mode) for the merged normal
 */
process GATK_MUTECT2_MITO_NORMAL {

    tag "$patient"
    cpus 18
    memory '32GB'
    time '45m'  // 2h is usually sufficient
    executor 'slurm'
    queue 'short'

    publishDir params.output_dir, mode: 'copy'

    input:
    tuple path(normal_bam), path(normal_bam_index)
    val patient
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)            // ref genome files
    path output_dir                                                 // location for output files
    tuple path(polymorphic_sites), path(polymorphic_sites_tbi)      // polymorphic sites
    tuple path(germline_resource), path(germline_resource_tbi)      // germline resource
    tuple path(panel_of_normals), path(panel_of_normals_idx)        // panel of normals

    output:
    tuple val("Nall"), path("${patient}_Nall_chrM.vcf.gz"), path("${patient}_Nall_chrM.vcf.gz.tbi"), path("${patient}_Nall_chrM.vcf.gz.stats"), path("${patient}_Nall_chrM.f1r2.tar.gz")

    script:
    """

    # run mutect2 for this region
    module load gcc/6.2.0 gatk/4.1.9.0

    gatk Mutect2 -R $ref_fasta \
        -L chrM \
        --mitochondria-mode true \
        -I ${normal_bam} \
        -O ${patient}_Nall_chrM.vcf.gz \
        --disable-read-filter MappingQualityAvailableReadFilter \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --f1r2-tar-gz ${patient}_Nall_chrM.f1r2.tar.gz \
        --native-pair-hmm-threads 16
    """
}


/*
 * filter mutect calls
 */
process FILTER_MITO_MUTECT_CALLS { 

    tag "$sample"
    cpus 8
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.mito_output_dir, mode: 'copy'

    input:
    val patient
    path tmp_dir
    tuple val(sample), path(raw_vcf), path(raw_vcf_tbi), path(raw_vcf_stats), path(raw_artifact_priors)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)  // ref genome files

    output:
    tuple path("${sample}_chrM_unfiltered_norm.vcf.gz"), path("${sample}_chrM_unfiltered_norm.vcf.gz.tbi"), path("${sample}_chrM_filtered.vcf.gz"), path("${sample}_chrM_filtered.vcf.gz.tbi")

    script:
    """
    module load gcc/6.2.0 gatk/4.1.9.0 bcftools/1.13

    # Learn read orientation model 
    gatk LearnReadOrientationModel -I ${raw_artifact_priors} -O ${sample}_read-orientation-model.tar.gz

    # FilterMutectCalls
    gatk FilterMutectCalls --mitochondria-mode true -R $ref_fasta -V $raw_vcf --orientation-bias-artifact-priors ${sample}_read-orientation-model.tar.gz -O ${sample}_chrM_unfiltered.vcf.gz

    # IndexFeatureFile
    gatk IndexFeatureFile -I ${sample}_chrM_unfiltered.vcf.gz --tmp-dir $tmp_dir

    # normalize VCF to split multi-allelic sites
    bcftools norm --multiallelics -both --fasta-ref $ref_fasta ${sample}_chrM_unfiltered.vcf.gz | bcftools view -I -O z -o ${sample}_chrM_unfiltered_norm.vcf.gz -

    # indexing normalized VCF
    gatk IndexFeatureFile -I ${sample}_chrM_unfiltered_norm.vcf.gz --tmp-dir $tmp_dir

    # filtering mutations
    bcftools view -i "%FILTER='PASS'" ${sample}_chrM_unfiltered_norm.vcf.gz | bcftools view -I -O z -o ${sample}_chrM_filtered.vcf.gz -

    # indexing filtered VCF
    gatk IndexFeatureFile -I ${sample}_chrM_filtered.vcf.gz --tmp-dir $tmp_dir
    """
}


/*
 * filter mutect calls
 */
process FILTER_MITO_MUTECT_CALLS_NORMAL { 

    tag "$sample"
    cpus 8
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.mito_output_dir, mode: 'copy'

    input:
    val patient
    path tmp_dir
    tuple val(sample), path(raw_vcf), path(raw_vcf_tbi), path(raw_vcf_stats), path(raw_artifact_priors)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)  // ref genome files

    output:
    tuple path("${sample}_chrM_unfiltered_norm.vcf.gz"), path("${sample}_chrM_unfiltered_norm.vcf.gz.tbi"), path("${sample}_chrM_filtered.vcf.gz"), path("${sample}_chrM_filtered.vcf.gz.tbi")

    script:
    """
    module load gcc/6.2.0 gatk/4.1.9.0 bcftools/1.13
 
    # Learn read orientation model 
    gatk LearnReadOrientationModel -I ${raw_artifact_priors} -O ${sample}_read-orientation-model.tar.gz

    # FilterMutectCalls
    gatk FilterMutectCalls --mitochondria-mode true -R $ref_fasta -V $raw_vcf --orientation-bias-artifact-priors ${sample}_read-orientation-model.tar.gz -O ${sample}_chrM_unfiltered.vcf.gz

    # IndexFeatureFile
    gatk IndexFeatureFile -I ${sample}_chrM_unfiltered.vcf.gz --tmp-dir $tmp_dir

    # normalize VCF to split multi-allelic sites
    bcftools norm --multiallelics -both --fasta-ref $ref_fasta ${sample}_chrM_unfiltered.vcf.gz | bcftools view -I -O z -o ${sample}_chrM_unfiltered_norm.vcf.gz -

    # indexing normalized VCF
    gatk IndexFeatureFile -I ${sample}_chrM_unfiltered_norm.vcf.gz --tmp-dir $tmp_dir

    # filtering mutations
    bcftools view -i "%FILTER='PASS'" ${sample}_chrM_unfiltered_norm.vcf.gz | bcftools view -I -O z -o ${sample}_chrM_filtered.vcf.gz -

    # indexing filtered VCF
    gatk IndexFeatureFile -I ${sample}_chrM_filtered.vcf.gz --tmp-dir $tmp_dir
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
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)  // ref genome files

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


process FILTER_RNA_EDIT_SITES {

    tag "$patient"
    cpus 8
    memory '16GB'
    time '10m'
    executor 'slurm'
    queue 'short'

    publishDir params.output_dir, mode: 'copy'

    input:
    val patient
    path filtered_vcf
    path filtered_vcf_tbi
    path rediportal

    output:
    path "${patient}_filtered_rediportal.vcf"

    script:
    """
    module load gcc/6.2.0 bcftools/1.13
    bcftools view -T ^${rediportal} ${filtered_vcf} -o ${patient}_filtered_rediportal.vcf
    
    """
}


process GENOTYPEBAM {
    tag "$sample"
    cpus 4
    memory '16GB'
    time '1h'
    executor 'slurm'
    queue 'short'

    publishDir params.maf_dir, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bam_index), path(barcode)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)
    path filtered_vcf

    output:
    tuple path("${sample}_ALT_matrix.txt"), path("${sample}_REF_matrix.txt"), path("${sample}_matrix_variants.txt")

    script:
    """
    /home/alg2264/nextflow/snrna_multisample/vartrix_linux --bam ${bam} --cell-barcodes ${barcode} --fasta ${ref_fasta} --vcf ${filtered_vcf} --scoring-method coverage --out-matrix ${sample}_ALT_matrix.txt --ref-matrix ${sample}_REF_matrix.txt --out-variants ${sample}_matrix_variants.txt --threads 4
    """
}


process GENOTYPEBAM_MITO {
    tag "$sample"
    cpus 4
    memory '16GB'
    time '10h'
    executor 'slurm'
    queue 'short'

    publishDir params.mito_output_dir, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bam_index), path(barcode)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)
    tuple path(filtered_vcf), path(filtered_vcf_csi)

    output:
    tuple path("${sample}_mtdna_ALT_matrix.txt"), path("${sample}_mtdna_REF_matrix.txt"), path("${sample}_mtdna_matrix_variants.txt")

    script:
    """
    /home/alg2264/nextflow/snrna_multisample/vartrix_linux --bam ${bam} --cell-barcodes ${barcode} --fasta ${ref_fasta} --vcf ${filtered_vcf} --scoring-method coverage --out-matrix ${sample}_mtdna_ALT_matrix.txt --ref-matrix ${sample}_mtdna_REF_matrix.txt --out-variants ${sample}_mtdna_matrix_variants.txt --threads 4
    """
}


process GENOTYPEBAM_MITO_NORMAL {
    tag "$sample"
    cpus 4
    memory '16GB'
    time '1h'
    executor 'slurm'
    queue 'short'

    publishDir params.mito_output_dir, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bam_index), path(barcode)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)
    tuple path(filtered_vcf), path(filtered_vcf_csi)
    path split_bam_dir

    output:
    tuple path("${sample}_normal_mtdna_ALT_matrix.txt"), path("${sample}_normal_mtdna_REF_matrix.txt"), path("${sample}_normal_mtdna_matrix_variants.txt")

    script:
    """
    /home/alg2264/nextflow/snrna_multisample/vartrix_linux --bam ${split_bam_dir}/${sample}_normal.bam --cell-barcodes ${split_bam_dir}/${sample}_normal_barcodes.txt_noheader.txt --fasta ${ref_fasta} --vcf ${filtered_vcf} --scoring-method coverage --out-matrix ${sample}_normal_mtdna_ALT_matrix.txt --ref-matrix ${sample}_normal_mtdna_REF_matrix.txt --out-variants ${sample}_normal_mtdna_matrix_variants.txt --threads 4
    """
}


process MERGE_MITO_VCFS {
    tag "$patient"
    cpus 4
    memory '8GB'
    time '5m'
    executor 'slurm'
    queue 'short'

    publishDir params.mito_output_dir, mode: 'copy'

    input:
    val patient
    path all_tumor_filtered_vcf
    path all_tumor_filtered_vcf_tbi
    path normal_filtered_vcf
    path normal_filtered_vcf_tbi

    output:
    tuple path("${patient}_filtered_vcf_merged.vcf.gz"), path("${patient}_filtered_vcf_merged.vcf.gz.csi")

    script:
    def all_tumor_filtered_vcf_line = all_tumor_filtered_vcf.collect { vcf -> "${vcf}" }.join(' ')

    """
    module load gcc/6.2.0 bcftools/1.13
    bcftools merge -m all -O z -o ${patient}_filtered_vcf_merged.vcf.gz ${all_tumor_filtered_vcf} ${normal_filtered_vcf}
    bcftools index ${patient}_filtered_vcf_merged.vcf.gz
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
    tuple val(sample), val(fq_prefix), val(sample_order)
    path filtered_vcf

    output:
    path "${sample}.maf"

    script:
    """
    module load gcc/9.2.0 bcftools/1.14 samtools/1.15.1

    # Replacement logic
    fixed_sample="$sample"
    case "\$fixed_sample" in
        "S12")
            fixed_sample="NaxID-ff-01-S12_snRNAseq_001"
            ;;
        "S14")
            fixed_sample="NaxID-ff-01-S14_snRNAseq_001"
            ;;
        "S9")
            fixed_sample="NaxID-ff-01-S9_snRNAseq_001"
            ;;
    esac

    # subset the multi-sample VCF for this sample
    bcftools view $filtered_vcf -s \$fixed_sample > ${sample}.vcf

    # run vcf2maf with the corrected tumor ID
    conda run -n vep perl /home/alg2264/repos/vcf2maf/vcf2maf.pl \\
        --input-vcf ${sample}.vcf \\
        --output-maf ${sample}.maf \\
        --tumor-id \$fixed_sample \\
        --remap-chain /home/alg2264/repos/vcf2maf/data/hg38_to_GRCh38.chain
    """
}



/*
 * Workflow
 */
workflow {

    // reference genome inputs
    ref_fasta = file(params.ref_fasta)
    ref_fai = file(params.ref_fai)
    ref_dict = file(params.ref_dict)
    ref_files = tuple(ref_fasta, ref_fai, ref_dict)

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

    // channel of genomic chunks
    genome_chunk_ch = Channel.fromPath(params.genome_chunks)
                    .splitCsv(header: false, sep: ",", strip: true)
                    .map { row -> tuple(row[0], row[1], row[2], row[3]) }


    // =============================
    // variant calling, filtering
    // =============================

    // run process to make all the expected directories for this patient
    make_dirs_ch = MAKE_DIRS(params.output_dir, params.maf_dir, params.mtbam_dir, params.tmp_dir, params.bam_dir)

    // split N reads (required for RNAseq)
    splitntrim_output = GATK_SPLITNTRIM(sample_ch, params.patient, params.bam_dir, params.split_bam_dir, ref_files)

    // Extracting each element into separate channels
    split_sample_ch = splitntrim_output.map { it[0] }
    split_tumor_bam_ch = splitntrim_output.map { it[1] }
    split_tumor_bam_index_ch = splitntrim_output.map { it[2] }
    split_normal_bam_ch = splitntrim_output.map { it[3] }
    split_normal_bam_index_ch = splitntrim_output.map { it[4] }

    // collect the bam files so that we can do multi-sample variant calling
    all_split_sample_ch = split_sample_ch.collect()
    all_split_tumor_bam_ch = split_tumor_bam_ch.collect()
    all_split_tumor_bam_index_ch = split_tumor_bam_index_ch.collect()
    all_split_normal_bam_ch = split_normal_bam_ch.collect()
    all_split_normal_bam_index_ch = split_normal_bam_index_ch.collect()

    // call mutations in multi-sample paired T/N mode
    mergenormal_output = GATK_MERGE_NORMAL(all_split_normal_bam_ch, all_split_normal_bam_index_ch, params.patient)
    
    // call mutations in multi-sample paired T/N mode
    mutect_output = GATK_MUTECT2(genome_chunk_ch, params.targets_bed, all_split_tumor_bam_ch, all_split_tumor_bam_index_ch, mergenormal_output, params.patient, params.output_dir, polymorphic_sites_files, germline_resource_files, panel_of_normals_files, ref_files)

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

    // remove RNA-editing sites
    filtered_rnaedit_output = FILTER_RNA_EDIT_SITES(params.patient, filtered_vcf_ch, filtered_vcf_tbi_ch, params.rediportal)

    // genotype each cell (one sample at a time) for the filtered variants
    genotypebam_output = GENOTYPEBAM(sample_ch, ref_files, filtered_rnaedit_output)

    // VCF2MAF
    vcf2maf_output = VCF2MAF(sample_ch, filtered_rnaedit_output)

    // MUTECT2_MITO_MODE
    mutect2_mtdna = GATK_MUTECT2_MITO(splitntrim_output, params.patient, ref_files, params.output_dir, polymorphic_sites_files, germline_resource_files, panel_of_normals_files)
    mutect2_mtdna_normal = GATK_MUTECT2_MITO_NORMAL(mergenormal_output, params.patient, ref_files, params.output_dir, polymorphic_sites_files, germline_resource_files, panel_of_normals_files)
    filter_mtdna_calls_output = FILTER_MITO_MUTECT_CALLS(params.patient, params.tmp_dir, mutect2_mtdna, ref_files)
    filter_mtdna_normal_calls_output = FILTER_MITO_MUTECT_CALLS_NORMAL(params.patient, params.tmp_dir, mutect2_mtdna_normal, ref_files)

    // collect all the tumor mtdna .vcfs and .vcf.tbi's, also the single normal vcf/vcf.tbi, then take their union
    tumor_mtdna_vcf_ch = filter_mtdna_calls_output.map { it[2] }
    tumor_mtdna_vcf_tbi_ch = filter_mtdna_calls_output.map { it[3] }
    all_tumor_mtdna_vcf = tumor_mtdna_vcf_ch.collect()
    all_tumor_mtdna_vcf_tbi = tumor_mtdna_vcf_tbi_ch.collect()
    normal_mtdna_vcf_ch = filter_mtdna_normal_calls_output.map { it[2] }
    normal_mtdna_vcf_tbi_ch = filter_mtdna_normal_calls_output.map { it[3] }
    merge_mito_vcfs_output = MERGE_MITO_VCFS(params.patient, all_tumor_mtdna_vcf, all_tumor_mtdna_vcf_tbi, normal_mtdna_vcf_ch, normal_mtdna_vcf_tbi_ch)

    // genotype each cell (one sample at a time) for the filtered variants
    mtdna_genotypebam_output = GENOTYPEBAM_MITO(sample_ch, ref_files, merge_mito_vcfs_output)
    mtdna_normal_genotypebam_output = GENOTYPEBAM_MITO_NORMAL(sample_ch, ref_files, merge_mito_vcfs_output, params.split_bam_dir)


}



