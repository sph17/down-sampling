version 1.0

import "Structs.wdl"

workflow collectOriginalReadCounts {

  #################################################################################
  ####        Required basic arguments for collectOriginalReadCounts              #
  #################################################################################
    
  input {
    File reference_fasta
    String downsample_docker
    File cram_or_bam_file
    File ref_fai
    File ref_dict

    String sample_ID
    File intervals_genome

    File? gatk4_jar_override
    String gatk_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_collect_counts

    Boolean run_count_coverage = true

  }

  parameter_meta {
    cram_or_bam_file: "CRAM/BAM file of interest"
    reference_fasta: ".fasta file with reference used to align bam or cram file"
  }

  meta {
    author: "Stephanie Hao"
    email: "shao@broadinstitute.org"
  }
    


  #################################################################################
  ####        Collect Read Counts with gatk                                       #
  #################################################################################
  if (run_count_coverage){
    call collectCountsCram {
    input : 
      cram = cram_or_bam_file,
      crai = cram_or_bam_file + ".crai",
      sample_ID = sample_ID,
      hg38_reference = reference_fasta,
      hg38_reference_fai = ref_fai,
      hg38_reference_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
    }
  }
  
  output {
    File? orginal_counts = collectCountsCram.counts_reads
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
    String counts_reads_filename_zip = "${counts_reads_filename}" + ".gz"
    
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

    gzip -c ~{counts_reads_filename} > ~{counts_reads_filename_zip}

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
    File counts_reads = counts_reads_filename_zip
  }
}

