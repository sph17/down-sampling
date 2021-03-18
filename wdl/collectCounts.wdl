version 1.0

import "Structs.wdl"

workflow collectReadCounts {

  #################################################################################
  ####        Required basic arguments for downsampling part 2                    #
  #################################################################################
    
  input {
    File reference_fasta
    File cram_sorted_file_2x
    File crai_file_2x
    File cram_sorted_file_4x
    File crai_file_4x
    File cram_sorted_file_6x
    File crai_file_6x
    File cram_sorted_file_8x
    File crai_file_8x
    
    File ref_fai
    File ref_dict

    String sample_ID
    File intervals_genome

    File? gatk4_jar_override
    String gatk_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_collect_counts

  }

  parameter_meta {
    cram_sorted_file: "The downsampled cram files"
    reference_fasta: ".fasta file with reference used to align bam or cram file"
  }

  meta {
    author: "Stephanie Hao"
    email: "shao@broadinstitute.org"
  }
    
  #################################################################################
  ####        Collect Read Counts with gatk                                       #
  #################################################################################
  call collectCountsCram as cc2x {
    input :
      intervals_genome = intervals_genome,
      cram = cram_sorted_file_2x,
      crai = crai_file_2x,
      sample_ID = sample_ID,
      hg38_reference = reference_fasta,
      hg38_reference_fai = ref_fai,
      hg38_reference_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  call collectCountsCram as cc4x {
    input :
      intervals_genome = intervals_genome,
      cram = cram_sorted_file_4x,
      crai = crai_file_4x,
      sample_ID = sample_ID,
      hg38_reference = reference_fasta,
      hg38_reference_fai = ref_fai,
      hg38_reference_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  call collectCountsCram as cc6x {
    input :
      intervals_genome = intervals_genome,
      cram = cram_sorted_file_6x,
      crai = crai_file_6x,
      sample_ID = sample_ID,
      hg38_reference = reference_fasta,
      hg38_reference_fai = ref_fai,
      hg38_reference_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  call collectCountsCram as cc8x {
    input :
      intervals_genome = intervals_genome,
      cram = cram_sorted_file_8x,
      crai = crai_file_8x,
      sample_ID = sample_ID,
      hg38_reference = reference_fasta,
      hg38_reference_fai = ref_fai,
      hg38_reference_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  output {
    File read_count_2x = cc2x.counts_reads
    File read_count_4x = cc4x.counts_reads
    File read_count_6x = cc6x.counts_reads
    File read_count_8x = cc8x.counts_reads
  }
}

task collectCountsCram {
  input {
    File intervals_genome
    File cram
    File crai
    File hg38_reference
    File hg38_reference_fai
    File hg38_reference_dict
    String sample_ID
    String gatk_docker
    File? gatk4_jar_override 
    RuntimeAttr? runtime_attr_override
  }

  # Runtime parameters adapted from gatk-sv "CollectCoverage.wdl"
  Int num_cpu = 1
  Int mem_size_gb = 6
  Int vm_disk_size = 50

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_overhead_gb = 2.0
  Int command_mem_mb = floor((mem_size_gb - mem_overhead_gb) * 1024)


  String base_filename = basename(cram, ".cram")
  String counts_reads_filename = "${sample_ID}.counts.tsv"
    
  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR="/root/gatk.jar"
    #Collects read counts and output as a counts.tsv file for downstream analyses
    gatk --java-options "-Xmx~{command_mem_mb}m" CollectReadCounts \
      -I ~{cram} \
      --read-index ~{crai} \
      -L ~{intervals_genome} \
      --interval-merging-rule OVERLAPPING_ONLY \
      --reference ~{hg38_reference} \
      --format TSV \
      -O ~{counts_reads_filename}

  >>>

  runtime {
    docker: gatk_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File counts_reads = counts_reads_filename
  }
}


