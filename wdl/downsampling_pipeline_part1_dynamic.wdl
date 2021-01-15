version 1.0

import "Structs.wdl"

workflow downSampling_01 {

  #################################################################################
  ####        Required basic arguments for downsampling part 1                    #
  #################################################################################
    
  input {
    File bam_or_cram_file
    File reference_fasta
    String downsample_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_bam_to_fastq    
  }

  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
    bam_file: ".bam file to search for SVs. bams are preferable, crams will be converted to bams."
    reference_fasta: ".fasta file with reference used to align bam or cram file"
  }

  meta {
    author: "Stephanie Hao"
    email: "shao@broadinstitute.org"
  }
    
  #################################################################################
  ####        Convert cram to bam                                                 #
  #################################################################################
  
  Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)

  if (!is_bam_) {
    call cramToBam {
      input :
        cram_file = bam_or_cram_file,
        reference_fasta = reference_fasta,
        downsample_docker = downsample_docker,
        runtime_attr_override = runtime_attr_cram_to_bam
    }
  }


  #################################################################################
  ####        Convert bam to fq1 & fq2                                            #
  #################################################################################
  call bamToFq { 
    input :
      bam_file = select_first([cramToBam.bam_file, bam_or_cram_file]),
      downsample_docker = downsample_docker,
      runtime_attr_override = runtime_attr_bam_to_fastq
  }

  output {
    File fastq_1 = bamToFq.fastq_file_1
    File fastq_2 = bamToFq.fastq_file_2
    File read_groups = cramToBam.read_groups_file
  }

}


task cramToBam {
    
  input {
    File cram_file
    File reference_fasta
    String downsample_docker
    RuntimeAttr? runtime_attr_override
  }

  File reference_index_file = reference_fasta + ".fai"
    
  String bam_file_name = basename(cram_file, ".cram") + ".bam"
  String read_groups_name = basename(cram_file, ".cram") + "_read_groups.txt"

  Int num_cpu = 4
  Float mem_size_gb = num_cpu * 4.0
  
  Float cram_inflate_ratio = 3.5
  Float disk_overhead = 20.0
  Float cram_size = size(cram_file, "GiB")
  Float bam_size = cram_inflate_ratio * cram_size
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(bam_size + ref_size + ref_index_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bam_file = bam_file_name
    File read_groups_file = read_groups_name
  }

  command {
    
    set -euo pipefail

    #extracts read groups for part 2 of pipeline
    samtools view -H ~{cram_file} > ~{read_groups_name}
    
    #converts cram files to bam files
    samtools view \
            -b \
            -h \
            -@ ~{num_cpu} \
            -T "~{reference_fasta}" \
            -o "~{bam_file_name}" \
            "~{cram_file}"
    }

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task bamToFq {
  
  input {
    File bam_file
    String downsample_docker
    RuntimeAttr? runtime_attr_override
  }

  String fastq_file_1_name = basename(bam_file, ".bam") + "_1.fastq"
  String fastq_file_2_name = basename(bam_file, ".bam") + "_2.fastq"

  Int num_cpu = 5
  Int mem_size_gb = 14
  Int vm_disk_size = 300

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }  

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File fastq_file_1 = fastq_file_1_name
    File fastq_file_2 = fastq_file_2_name
  }

  command <<<
    set -euo pipefail

    #converts bam file to paired fastq files
    picard SamToFastq \
            I=~{bam_file} \
            FASTQ=~{fastq_file_1_name} \
            SECOND_END_FASTQ=~{fastq_file_2_name}
  >>>

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}



