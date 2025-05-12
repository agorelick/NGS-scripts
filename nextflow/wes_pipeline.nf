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


// Create tuple channel: (sample, fq_prefix, order)
def sample_ch = Channel.from( rows.collect { [it.sample, it.fq_prefix, it.order] } )

// Set scalar params using first row with non-null value
def first_full_row = rows.find { it.patient }  // or any other required field
params.patient        = first_full_row.patient
params.normal_sample = first_full_row.normal_sample
params.sex            = first_full_row.sex
params.fq_dir         = first_full_row.fq_dir
params.nextflow_dir    = first_full_row.nextflow_dir
params.output_dir   = first_full_row.output_dir

// dynamically generated parameters
params.bam_dir = "${params.output_dir}/${params.patient}/bams"
params.maf_dir = "${params.output_dir}/${params.patient}/mafs"
params.mosdepth_dir = "${params.output_dir}/${params.patient}/mosdepth"
params.mtbam_dir = "${params.output_dir}/${params.patient}/mtbams"
params.haplocheck_dir = "${params.output_dir}/${params.patient}/haplocheck"
params.ascat_dir = "${params.output_dir}/${params.patient}/ascat"
params.tmp_dir = "${params.nextflow_dir}/${params.patient}/tmp_files"


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
params.genome_chunks = "/n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38_50Mbchunks.csv"

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
    val bam_dir
    val maf_dir
    val mosdepth_dir
    val mtbam_dir
    val haplocheck_dir
    val ascat_dir
    val tmp_dir

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

    touch mkdir_done.txt
    """
}


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
    tuple val(sample), val(fq_prefix), val(sample_order)
    val fq_dir
    val make_dirs_ch

    output:
    tuple val(sample), val(sample_order), path("${sample}_trimmed_R1.fastq.gz"), path("${sample}_trimmed_R2.fastq.gz")

    script:
    """
    module load gcc/9.2.0 python/3.8.12 cutadapt/4.1

    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --cores=8 \\
        -o ${sample}_trimmed_R1.fastq.gz -p ${sample}_trimmed_R2.fastq.gz \\
        ${fq_dir}/${fq_prefix}_R1_001.fastq.gz ${fq_dir}/${fq_prefix}_R2_001.fastq.gz
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
    tuple val(sample), val(sample_order), path(trimmed_fq1), path(trimmed_fq2)
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_raw.sam")

    script:
    """
    module load gcc/6.2.0 bwa/0.7.15

    bwa mem -M -t 8 -R '@RG\\tID:${sample_order}\\tSM:${sample}\\tPL:Illumina' ${ref_fasta} ${trimmed_fq1} ${trimmed_fq2} > ${sample}_raw.sam

    """
}


/*
 * Run picard SortSam
 */
process PICARD_SORTSAM {

    tag "$sample"
    cpus 8
    memory '72GB'
    time '3h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), path(raw_sam)
    path bam_dir
    path tmp_dir
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_sorted.bam"), path("${sample}_sorted.bai")

    script:
    """
    mkdir -p ${bam_dir}

    module load picard/2.27.5
    java -jar $PICARD/picard.jar SortSam --INPUT ${raw_sam} --OUTPUT ${sample}_sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true --TMP_DIR ${tmp_dir} -R ${ref_fasta}  
    """
}



/*
 * Run GATK MarkDuplicates
 */
process GATK_MARKDUP {

    tag "$sample"
    cpus 8
    memory '72GB'
    time '12h'
    executor 'slurm'
    queue 'short'

    input:
    tuple val(sample), path(sorted_bam), path(sorted_bai)
    path bam_dir
    path tmp_dir
    tuple path(ref_fasta), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_fai), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample), path("${sample}_marked_dup.bam"), path("${sample}_marked_dup.bai")

    script:
    """

    mkdir -p ${bam_dir}
    module load picard/2.27.5
    java -Xmx60g -jar $PICARD/picard.jar MarkDuplicates --INPUT ${sorted_bam} --OUTPUT ${sample}_marked_dup.bam --ASSUME_SORT_ORDER coordinate --CREATE_INDEX true --TMP_DIR ${tmp_dir} -R ${ref_fasta} --METRICS_FILE ${bam_dir}/${sample}_marked_dup_metrics.txt

 
    """
}

/*
 * Run GATK BaseRecalibrator
 */
process GATK_BASERECAL {

    tag "$sample"
    cpus 16
    memory '24GB'
    time '12h'
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
    memory '32GB'
    time '3h'
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
    echo "foo"
    gatk ApplyBQSR -R ${ref_fasta} -I ${markdup_bam} --bqsr-recal-file ${recal_data_table} -O ${sample}.bam --tmp-dir ${tmp_dir}

    """
}


/*
 * Run GATK Mutect2 for multi-sample tumor/normal variant calling
 */
process GATK_MUTECT2 {

    tag "$region"
    cpus 18
    memory '32GB'
    time '12h'  // 2h is usually sufficient
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
    tuple val(sample), val(fq_prefix), val(sample_order)
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

    ## run mosdepth for this sample's bam file
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
    module load gcc/9.2.0 bcftools/1.14 samtools/1.15.1
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

    module unload R/4.3.1
    conda run -n CNalign Rscript /home/alg2264/repos/CNalign/scripts/run_ascat_prepareHTS.R \
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
    path "${patient}_CNalign_obj.rds"
    path "${patient}_CNalign_obj_mpcf.rds"
    path "${patient}_CNalign_obj_mpcf_hisens.rds"



    script:
    """
    echo "Generating CNalign data object ..."
    conda run -n CNalign /home/alg2264/miniconda3/envs/CNalign/bin/Rscript /home/alg2264/repos/CNalign/scripts/merge_alleleCounter_data.R \
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

    // channel of genomic chunks
    genome_chunk_ch = Channel.fromPath(params.genome_chunks)
                    .splitCsv(header: false, sep: ",", strip: true)
                    .map { row -> tuple(row[0], row[1], row[2], row[3]) }


    // =============================
    // bam preprocessing
    // =============================

    // run process to make all the expected directories for this patient
    make_dirs_ch = MAKE_DIRS(params.bam_dir, params.maf_dir, params.mosdepth_dir, params.mtbam_dir, params.haplocheck_dir, params.ascat_dir, params.tmp_dir)

    // Adapter trimming
    trim_adapt_output = TRIM_ADAPTERS(sample_ch, params.fq_dir, make_dirs_ch)

    // BWA-MEM alignment
    bwa_mem_output = BWA_MEM(trim_adapt_output, ref_files)

    // Run SortSam on valid samples
    sortsam_output = PICARD_SORTSAM(bwa_mem_output, params.bam_dir, params.tmp_dir, ref_files)

    // Run MarkDuplicates on sorted bam files
    markdup_output = GATK_MARKDUP(sortsam_output, params.bam_dir, params.tmp_dir, ref_files)

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

    // mosdepth (QC)
    mosdepth_output = MOSDEPTH(applybqsr_output, params.targets_bed)

    // slice mtDNA
    slice_mtdna_output = SLICE_MTDNA(applybqsr_output, params.mt_label)

    // Extracting each element into separate channels
    mtsample_ch = slice_mtdna_output.map { it[0] }
    mtbam_ch = slice_mtdna_output.map { it[1] }
    mtbam_index_ch = slice_mtdna_output.map { it[2] }
    all_mtbam_ch = mtbam_ch.collect()
    all_mtbam_index_ch = mtbam_index_ch.collect()

    // 1: Separate tumor and normal
    tumor_input_ch = applybqsr_output.filter { sample, bam, bai -> sample != params.normal_sample }
    normal_input_ch = applybqsr_output.filter { sample, bam, bai -> sample == params.normal_sample }

    // 2. Flatten into all pairwise combinations
    prep_input_ch = tumor_input_ch.combine(normal_input_ch)

    // run for each tuple (which is a combination of one tumor and the same repeated normal)
    prep_data_output = PREP_CNA_DATA(prep_input_ch, params.patient, params.sex, params.build, params.targets_bed, params.allelecounter_exe, params.alleles_prefix, params.loci_prefix)
    all_allelecounter_files_ch = prep_data_output.collect()
    cnalign_output = GET_CNALIGN_OBJ(all_allelecounter_files_ch, params.normal_sample, params.patient, params.sex, params.build, params.gccontentfile, params.replictimingfile)
}



