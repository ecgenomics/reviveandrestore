#!/usr/bin/env nextflow

/*
  _____            _                      _____           _                 
 |  __ \          (_)             ___    |  __ \         | |                
 | |__) |_____   _____   _____   ( _ )   | |__) |___  ___| |_ ___  _ __ ___ 
 |  _  // _ \ \ / / \ \ / / _ \  / _ \/\ |  _  // _ \/ __| __/ _ \| '__/ _ \
 | | \ \  __/\ V /| |\ V /  __/ | (_>  < | | \ \  __/\__ \ || (_) | | |  __/
 |_|  \_\___| \_/ |_| \_/ \___|  \___/\/ |_|  \_\___||___/\__\___/|_|  \___|
                                                                            
                                                                            
  Module name: main.nf
                                                                           
  Endogenous or hDNA Quantification to prepare equi-endogenous pools for capture 
  Based on https://github.com/claudefa/Mapping_Capture_Experiments
*/


// Syntax

nextflow.enable.dsl=2

// Inputs
///////////////////////////////////////////////////////

params.reads = "input/test_R{1,2}.fastq"
params.genome = "input/hg19.fa.masked"
params.bedfile = "input/tuf.bed"
// Folders
///////////////////////////////////////////////////////

//    Session
params.session = baseDir + "/output"

// Temp
params.temp = params.session + "/temp"

//    Demultiplexed fastq
params.dm = params.session + "/fastqs/demultiplex/"

//    Trimmed adapters outputs
params.tr_out = params.session + "/fastqs/trimmed"                        // Output trimmed out 

//    Mapping HG19 hDNA
params.bam_hg19 = params.session + "/bam_hg19"

//    Remove duplicates
params.bam_hg19_no_dups = params.session + "/bam_hg19_no_dups"

//    Quality filtering
params.BAM_hg19_filtered = params.session + "/BAM_hg19_filtered"

//    On target
params.BAM_ontarget = params.session + "/BAM_ontarget"


// Software (warning: this needs to change with the singularity)
///////////////////////////////////////////////////////

//    sabre
params.sabre = baseDir + "/software/sabre/sabre"

//    picardtools
params.picard = baseDir + "/software/picard/build/libs/picard.jar"

//    bwa
params.bwa = baseDir + "/software/bwa/bwa"

//    samtools
params.samtools = baseDir + "/software/samtools-1.20/samtools"

//    bedtools
params.bedtools = baseDir + "/software/bedtools2/bin/"

// Step 0: Multiplexing

process multiplexing {
    tag "$sample_id"
    publishDir "${params.dm}", mode: "copy", overwrite : true

    input:
    tuple val(sample_id), path(f), path(r)

    output:
    tuple val(sample_id), path("${sample_id}_demultiplexed_R1.fastq"), path("${sample_id}_demultiplexed_R2.fastq")

    script:
    """
    ${params.sabre} pe -m 2 -c -f ${f} -r ${r} -b  ${baseDir}/src/barcode_files/Pool_Example1_barcode_file.txt \
    -u  ${sample_id}_demultiplexed_R1.fastq -w  ${sample_id}_demultiplexed_R2.fastq
  
    """
}

// Step 1: Trimming

process trim{

    publishDir "${params.tr_out}", mode: "copy", overwrite : true

    input:
    tuple val(sample_id), path(f), path(r)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq"), path("${sample_id}_trimmed_R2.fastq"), path("${sample_id}_trimmed_unpaired_R1.fastq"), path("${sample_id}_trimmed_unpaired_R2.fastq"), path("${sample_id}.merged.fastq"), path("${sample_id}.report.html")

    script:
    """
    fastp -i ${f} -I ${r} \
    --merge --merged_out ${sample_id}.merged.fastq \
    --out1 ${sample_id}_trimmed_R1.fastq --out2 ${sample_id}_trimmed_R2.fastq \
    --unpaired1 ${sample_id}_trimmed_unpaired_R1.fastq --unpaired2 ${sample_id}_trimmed_unpaired_R2.fastq \
    --overlap_len_require 11 --detect_adapter_for_pe -p \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --trim_poly_g --poly_g_min_len 10 --length_required 30 \
    --html ${sample_id}.report.html

    """

}

// Step 2: [1] Single end alignment

process align_se{

    publishDir "${params.bam_hg19}", mode: "copy", overwrite : true

    input:
    tuple val(sample_id), path(trm1), path(trm2), path(unp1), path(unp2), path(mrgd), path(html)

    output:
    tuple val(sample_id), path("${sample_id}_se.bam")
    script:
    """
    ${params.bwa} mem -t 4 -M ${params.genome} ${mrgd} |
    ${params.samtools} view -Sbh - |
    java -Djava.io.tmpdir=${params.temp} \
    -jar ${params.picard} AddOrReplaceReadGroups \
    I=/dev/stdin O=${sample_id}_se.bam \
    SORT_ORDER=coordinate QUIET=TRUE COMPRESSION_LEVEL=9 \
    MAX_RECORDS_IN_RAM=150000 RGLB=${sample_id} \
    RGPL=Illumina RGPU=${sample_id} RGSM=${sample_id} \
    RGID=${sample_id} VALIDATION_STRINGENCY=SILENT
  
    """
}

// Step 2: [2] Paired ends alignment

process align_pe{

    publishDir "${params.bam_hg19}", mode: "copy", overwrite : true

    input:
    tuple val(sample_id), path(f), path(r), path(unp1), path(unp2), path(mrgd), path(html)

    output:
    tuple val(sample_id), path("${sample_id}_pe.bam")

    script:
    """
    ${params.bwa} mem -t 4 -M ${params.genome} ${f} ${r} | 
    awk 'substr(\$0,1,1)==\"@\" || (\$9>= 35 && \$9<=500) || (\$9<=-35 && \$9>=-500)' |
    ${params.samtools} view -Shb - |
    java -Djava.io.tmpdir=${params.temp} -jar ${params.picard} AddOrReplaceReadGroups I=/dev/stdin O=${sample_id}_pe.bam \
    SORT_ORDER=coordinate QUIET=TRUE COMPRESSION_LEVEL=9 MAX_RECORDS_IN_RAM=150000 \
    RGLB=${sample_id} RGPL=Illumina RGPU=${sample_id} RGSM=${sample_id} RGID=${sample_id}\
    VALIDATION_STRINGENCY=SILENT
  
    """
}

// Step 3 : [1] Merge single- and paired- end

process merge_sepe{

    publishDir "${params.bam_hg19}", mode: "copy", overwrite : true

    input:
    tuple val(sample_id), path(s), path(p)

    output:
    tuple val(sample_id), path("${sample_id}_srtd.bam")

    script:
    """
    java -Xmx8g -Djava.io.tmpdir=${params.temp} -jar ${params.picard} MergeSamFiles I=${p} I=${s} O=${sample_id}_mrgd.bam; \
	  ${params.samtools} sort -T ${sample_id}.mg -o ${sample_id}_srtd.bam ${sample_id}_mrgd.bam; \
	  ${params.samtools} index ${sample_id}_srtd.bam;
  
    """
}

// Step 3 : [2] remove the duplicates

process rmdup{

    publishDir "${params.bam_hg19_no_dups}", mode: "copy", overwrite : true

    input:
    tuple val(sample_id), path(bamsorted)

    output:
    tuple val(sample_id), path("${sample_id}_rmdups.bam"), path("${sample_id}_rmdups.bam.stats")

    script:
    """
    java -Xmx8g -Djava.io.tmpdir=${params.temp} -jar ${params.picard} MarkDuplicates \
    I=${bamsorted} \
    O=${sample_id}_rmdups.bam \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 METRICS_FILE=${sample_id}_rmdups.bam.stats  \
    REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_RECORDS_IN_RAM=1500000 \
    VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true COMPRESSION_LEVEL=9
    """
}

// Step 4 : filter the quality

process filterQuality {

  publishDir "${params.BAM_hg19_filtered}", mode: "copy", overwrite : true

  input:
  tuple val(sample_id), path(bs_nodups), path(bs_stats)

  output:
  tuple val(sample_id), path("${sample_id}_rmdups.sorted.bam"), path("${sample_id}_rmdups.qual.bam")
  script:
  """
  ${params.samtools} sort -o ${sample_id}_rmdups.sorted.bam ${bs_nodups}; \
	${params.samtools} index ${bs_nodups};
	${params.samtools} view -h -q 30 -F 256 ${sample_id}_rmdups.sorted.bam | samtools view -hb > ${sample_id}_rmdups.qual.bam; \
	${params.samtools} index ${sample_id}_rmdups.qual.bam;
  """
}

// Step 5 : On target Hg19

process ontarget {

  publishDir "${params.BAM_ontarget}", mode: "copy", overwrite : true

  input:
  tuple val(sample_id), path(nodups_sorted), path(nodups_qual)
  
  output:
  tuple val(sample_id), path(nodups_sorted)

  script:
  """
  ${params.bedtools}/intersectBed -abam  ${nodups_qual} -b ${params.bedfile} > ${sample_id}_ontarget.bam
  """
}

// Step 6 : [1] Stats


// THE WORKFLOW
///////////////////////////////////////////////////////

workflow {

    // Inputs
    reads_ch = Channel.fromFilePairs(params.reads, flat: true)

    // Process launch
    demulti = multiplexing(reads_ch)              // Step 0 : demultiplexing
    trimmed = trim(demulti)                       // Step 1 : trimming
    aligned_se = align_se(trimmed)                // Step 2 : [1] Pair end alignment
    aligned_pe = align_pe(trimmed)                // Step 2 : [2] Paired ends alignment

    allbams = aligned_se.join(aligned_pe)         
    merged = merge_sepe(allbams)                  // Step 3 : [1] Merge single- and paired- end
    nodups = rmdup(merged)                        // Step 3 : [2] remove the duplicates
    filtered = filterQuality(nodups)              // Step 4 : filter
    //ontrg = ontarget(nodups)                     // Step 5 : On target
}  

