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

// Create tuple channel: (sample, offtarget_bam_file, order)
def sample_ch = Channel.from( rows.collect { [it.sample, it.offtarget_bam_file, it.offtarget_bai_file, it.order] } )

// Set scalar params using first row with non-null value
def first_full_row = rows.find { it.patient }  // or any other required field
params.patient        = first_full_row.patient
params.normal_sample = first_full_row.normal_sample
params.sex            = first_full_row.sex
params.nextflow_dir    = first_full_row.nextflow_dir
params.output_dir   = first_full_row.output_dir

// dynamically generated parameters
params.offtarget_fq_dir = "${params.output_dir}/${params.patient}/fq_offtarget"
params.bam_dir = "${params.output_dir}/${params.patient}/bams"
params.mosdepth_dir = "${params.output_dir}/${params.patient}/mosdepth"
params.mtbam_dir = "${params.output_dir}/${params.patient}/mtbams"
params.haplocheck_dir = "${params.output_dir}/${params.patient}/haplocheck"
params.ascat_dir = "${params.output_dir}/${params.patient}/ascat"
params.glimpse_dir = "${params.output_dir}/${params.patient}/glimpse"

/*
 * Additional pipeline parameters (use for all WES data)
 */

// Genome reference files
params.build = "hg38"
params.mt_label = "chrM"
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
params.targets_bed = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38.bed"

// references for GLIMPSE2
params.glimpse_refsplit_dir = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/GLIMPSE_GRCh38/reference_panel/split"
params.glimpse_chunks_dir = "/home/alg2264/repos/NGS-scripts/nextflow"

// ascat/CNAlign params
params.allelecounter_exe = "/home/alg2264/miniconda3/envs/CNalign/bin/alleleCounter"
params.alleles_prefix = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_allelesAll_hg38/G1000_alleles_hg38_chr"
params.loci_prefix = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_lociAll_hg38/G1000_loci_GRCh38_chr"
params.gccontentfile = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/GC_G1000_hg38.txt"
params.replictimingfile = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/RT_G1000_hg38.txt"

println "params.nextflow_dir: ${params.nextflow_dir}"
println "params.output_dir: ${params.output_dir}"
println "params.bam_dir: ${params.bam_dir}"



process MAKE_DIRS {
    tag "mkdirs"
    executor 'local'

    input:
    val offtarget_fq_dir
    val bam_dir
    val mosdepth_dir
    val mtbam_dir
    val haplocheck_dir
    val ascat_dir
    val glimpse_dir

    output:
    path "mkdir_done.txt"

    """
    mkdir -p ${offtarget_fq_dir}
    mkdir -p ${bam_dir}
    mkdir -p ${mosdepth_dir}
    mkdir -p ${mtbam_dir}
    mkdir -p ${haplocheck_dir}
    mkdir -p ${ascat_dir}
    mkdir -p ${glimpse_dir}


    touch mkdir_done.txt
    """
}


/*
 * Convert off-target .bam to a pair of fastq-files
 */
process CONVERT_BAM_TO_FQ_PAIR {

    tag "$sample"
    cpus 8
    memory '16GB'
    time '10m'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), val(offtarget_bam_file), val(offtarget_bai_file), val(sample_order)
    val make_dirs_ch

    output:
    tuple val(sample), val(sample_order), path("${sample}_offtarget_R1.fastq.gz"), path("${sample}_offtarget_R2.fastq.gz")

    publishDir params.offtarget_fq_dir, mode: 'copy'

    script:
    """
    module load gcc/14.2.0 samtools/1.21
    samtools fastq -@ 8 -1 ${sample}_offtarget_R1.fastq.gz -2 ${sample}_offtarget_R2.fastq.gz -0 /dev/null -s /dev/null -n ${offtarget_bam_file}
    """
}


/*
 * Trim Illumina Universal Adapters
 */
process TRIM_ADAPTERS {

    tag "$sample"
    cpus 8
    memory '16GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), val(sample_order), val(untrimmed_fq1), val(untrimmed_fq2)

    output:
    tuple val(sample), val(sample_order), path("${sample}_trimmed_R1.fastq.gz"), path("${sample}_trimmed_R2.fastq.gz")

    script:
    """
    conda run -n cutadapt cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --cores=8 \\
        -o ${sample}_trimmed_R1.fastq.gz -p ${sample}_trimmed_R2.fastq.gz \\
        ${untrimmed_fq1} ${untrimmed_fq2}
    """
}


/*
 * Run BWA-MEM
 */
process BWA_MEM {

    tag "$sample"
    cpus 8
    memory '16GB'
    time '1h'
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

    bwa mem -M -t 8 -R '@RG\\tID:${sample_order}\\tSM:${sample}\\tPL:Illumina' ${ref_fasta} ${trimmed_fq1} ${trimmed_fq2} > ${sample}_raw.sam

    """
}


/*
 * Run samtools sort
 */
process SAMTOOLS_SORT {

    tag "$sample"
    cpus 8
    memory '36GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.bam_dir, mode: 'copy'

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
    """
}



/*
 * Run GATK MarkDuplicates
 */
process GATK_MARKDUP {

    tag "$sample"
    cpus 8
    memory '64GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.bam_dir, mode: 'copy'

    input:
    tuple val(sample), path(sorted_bam), path(sorted_bai)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_marked_dup.bam"), path("${sample}_marked_dup.bam.bai"), path("${sample}_marked_dup_metrics.txt")

    script:
    """
    conda run -n gatk_4.6.1.0 gatk --java-options "-Xmx56g" MarkDuplicatesSpark \
        -I ${sorted_bam} \
        -O ${sample}_marked_dup.bam \
        -M ${sample}_marked_dup_metrics.txt \
        -R ${ref_fasta} \
        --create-output-bam-index true \
        --spark-master "local[8]"
    """
}



/*
 * Run GATK BaseRecalibrator
 */
process GATK_BASERECAL {

    tag "$sample"
    cpus 16
    memory '24GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.bam_dir, mode: 'copy'

    input:
    tuple val(sample), path(markdup_bam), path(markdup_bam_index), path(markdup_metrics)
    tuple path(polymorphic_sites), path(polymorphic_sites_tbi)    
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path(markdup_bam), path(markdup_bam_index), path("${sample}_recal_data.table")

    script:
    """
    conda run -n gatk_4.6.1.0 gatk BaseRecalibrator -I ${markdup_bam} -R ${ref_fasta} --known-sites ${polymorphic_sites} -O ${sample}_recal_data.table
    """
}

/*
 * Run GATK ApplyBQSR
 */
process GATK_APPLYBQSR {

    tag "$sample"
    cpus 16
    memory '32GB'
    time '1h'
    executor 'slurm'
    queue 'short'

    publishDir params.bam_dir, mode: 'copy'

    input:
    tuple val(sample), path(markdup_bam), path(markdup_bam_index), path(recal_data_table)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bai")

    script:
    """
    conda run -n gatk_4.6.1.0 gatk ApplyBQSR -R ${ref_fasta} -I ${markdup_bam} --bqsr-recal-file ${recal_data_table} -O ${sample}.bam

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
    samtools view -hb ${sample}.bam -F 3840 -q 20 ${mt_label} > ${sample}_mt.bam
    samtools index ${sample}_mt.bam
    """
}


/*
 * Run haplocheck
 */
process HAPLOCHECK {

    tag "$patient"
    cpus 2
    memory '8GB'
    time '2h'
    executor 'slurm'
    queue 'short'

    input:
    val patient
    path all_mtbams
    path all_mtbam_indices
    path haplocheck_dir

    output:
    path "${haplocheck_dir}/${patient}"

    script:
    """
    ## run haplocheck on all bams in the bams/ directory
    cloudgene run haplocheck@1.3.2 --files $PWD --format bam --output ${haplocheck_dir}/${patient} --threads 2
    """
}


process GLIMPSE2_PHASE {

    tag "${chr}"
    cpus 20
    memory '20GB'
    time '2h'
    executor 'slurm'
    queue 'short'

    publishDir params.glimpse_dir, mode: 'copy'

    input:
    //val chr
    //val normal_sample
    //path normal_bam
    //path normal_index
    tuple val(chr), val(normal_sample), path(normal_bam), path(normal_index)

    output:
    path "chr${chr}_chunks"  // directory with all .bcf for this chromosome

    script:
    """
    set -euo pipefail

    REF="${params.glimpse_refsplit_dir}/1000GP.chr${chr}"
    REGIONS_TXT="${params.glimpse_chunks_dir}/glimpse_chunks_per_chr/chr${chr}_regions.tsv"    
    PATIENT="${params.patient}"
    mkdir -p chr${chr}_chunks

    # Loop over the regions text passed from Nextflow
    while IFS=\$'\\t' read -r CHR START END CHRINT; do
        [ -z "\$START" ] && continue
        GLIMPSE2_phase_static \\
          --bam-file "${normal_bam}" \\
          --reference "\${REF}_chr${chr}_\${START}_\${END}.bin" \\
          --output "chr${chr}_chunks/\${PATIENT}_${normal_sample}_imputed_${chr}_\${START}_\${END}.bcf" \\
          --threads 20
    done < "\${REGIONS_TXT}"
    """
}


process GLIMPSE2_LIGATE {

    tag "${params.patient}"
    cpus 1
    memory '20GB'
    time '30m'
    executor 'slurm'
    queue 'short'

    publishDir params.glimpse_dir, mode: 'copy'

    input:
    path chr_chunk_dirs

    output:
    tuple path("${params.patient}_${params.normal_sample}_glimpse2_ligated.bcf"), path("${params.patient}_${params.normal_sample}_glimpse2_ligated_cleaned.bcf"), path("${params.patient}_${params.normal_sample}_glimpse2_ligated.bcf.csi"), path("${params.patient}_${params.normal_sample}_glimpse2_ligated_cleaned.bcf.csi")

    script:
    """ 

    module load gcc/14.2.0 bcftools/1.21
    
    ## these files will be generated by this script
    bcf_file="${params.patient}_${params.normal_sample}_glimpse2_ligated.bcf"
    bcf_cleaned_file="${params.patient}_${params.normal_sample}_glimpse2_ligated_cleaned.bcf"
    
    rm -f list_imputed_files.txt
    touch list_imputed_files.txt

    for d in ${chr_chunk_dirs}; do
        ## list the glimpse .bcf files for each chunk in genomic order to a text file
        ls -1v \${d}/*.bcf >> list_imputed_files.txt
    done

    ## Run GLIMPSE_ligate to combine the phased data for each genomic chunk
    GLIMPSE2_ligate_static --input list_imputed_files.txt --output \$bcf_file --threads 1

    ## use bcftools to subset the data for high-confidence heterozygous SNPs
    bcftools view -g het -i 'FORMAT/GP>=0.90' \$bcf_file | bcftools annotate -x FORMAT/DS,FORMAT/GP -O b - > \$bcf_cleaned_file
    bcftools index \$bcf_cleaned_file
    """
}




/*
 * Slice bam for chrM/MT for haplocheck
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
 * Get data object for CNalign (copied from WES, update this for bin-level)
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
    path "${patient}_CNalign_obj.rds"
    path "${patient}_CNalign_obj_mpcf.rds"
    path "${patient}_CNalign_obj_mpcf_hisens.rds"



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


    // =============================
    // bam preprocessing
    // =============================

    // run process to make all the expected directories for this patient
    make_dirs_ch = MAKE_DIRS(params.offtarget_fq_dir, params.bam_dir, params.mosdepth_dir, params.mtbam_dir, params.haplocheck_dir, params.ascat_dir, params.glimpse_dir)

    // convert offtarget bam file to a paired of FASTQ files
    bam_to_fq_ch = CONVERT_BAM_TO_FQ_PAIR(sample_ch, make_dirs_ch)

    // Adapter trimming
    trim_adapt_output = TRIM_ADAPTERS(bam_to_fq_ch)

    // BWA-MEM alignment
    bwa_mem_output = BWA_MEM(trim_adapt_output, ref_files)

    // Run SortSam on valid samples
    sortsam_output = SAMTOOLS_SORT(bwa_mem_output, ref_files)

    // Run MarkDuplicates on sorted bam files
    markdup_output = GATK_MARKDUP(sortsam_output, ref_files)

    // BaseRecalibrator
    baserecal_output = GATK_BASERECAL(markdup_output, polymorphic_sites_files, ref_files)

    // ApplyBQSR
    applybqsr_output = GATK_APPLYBQSR(baserecal_output, ref_files)

    // Extracting each element into separate channels
    preprocessed_sample_ch = applybqsr_output.map { it[0] }
    preprocessed_bam_ch = applybqsr_output.map { it[1] }
    preprocessed_bam_index_ch = applybqsr_output.map { it[2] }

    // Collect bams/indices
    all_bams_ch = preprocessed_bam_ch.collect()
    all_bam_indices_ch = preprocessed_bam_index_ch.collect()


    // =============================
    // run GLIMPSE2 on the normal bam 
    // =============================

    // 1) normal sample channel
    normal_input_ch = applybqsr_output.filter { sample, bam, bai -> sample == params.normal_sample }
    normal_sample_ch = normal_input_ch.map { it[0] }
    normal_bam_ch = normal_input_ch.map { it[1] }
    normal_bai_ch = normal_input_ch.map { it[2] }
 
    // 2) Chromosome channel
    chrom_ch = Channel.fromList( (1..22).collect { it.toString() } + 'X' )
    
    // 3) combine normal channel and chromosome channel
    normal_chr_ch = chrom_ch.combine(normal_input_ch)
    
    // run glimpse phase for each chromosome and collect the directories for their BCF files
    glimpse_phase_output = GLIMPSE2_PHASE(normal_chr_ch)
    grouped_glimpse_dirs_ch = glimpse_phase_output.collect()    
    glimpse_ligate_output = GLIMPSE2_LIGATE(grouped_glimpse_dirs_ch)

}



