version 1.0

import "Structs.wdl"

workflow compressCounts {

  #################################################################################
  ####        Required basic arguments for downsampling part 2                    #
  #################################################################################
    
  input {
    File counts_2x
    File counts_4x
    File counts_6x
    File counts_8x

    String sample_ID

    File? gatk4_jar_override
    String gatk_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_collect_counts

  }

  parameter_meta {
    counts_x: "collectCounts Files"
  }

  meta {
    author: "Stephanie Hao"
    email: "shao@broadinstitute.org"
  }
    
  #################################################################################
  ####        Collect Read Counts with gatk                                       #
  #################################################################################
  call compressCounts as cc2x {
    input :
      count_reads_filename = counts_2x,
      sample_ID = sample_ID,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  call compressCounts as cc4x {
    input :
      count_reads_filename = counts_6x,
      sample_ID = sample_ID,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  call compressCounts as cc6x {
    input :
      count_reads_filename = counts_6x,
      sample_ID = sample_ID,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  call compressCounts as cc8x {
    input :
      count_reads_filename = counts_8x,
      sample_ID = sample_ID,
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

task compressCounts {
  input {
    File count_reads_filename
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

  String count_reads_filename_zip = "${sample_ID}.counts.tsv.gz"
    
  command <<<
    set -euo pipefail

    gzip -c ~{count_reads_filename} > ~{count_reads_filename_zip}

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
    File counts_reads = count_reads_filename_zip
  }
}


