version 1.0

workflow downSampling_01 {

    #################################################################################
    ####        Required basic arguments for downsampling part 1                    #
    #################################################################################
    
    input {
        File bam_or_cram_file
        File reference_fasta
        Int? cram_to_bam_disk_size
        String? cram_to_bam_mem_size
        String downsample_docker

        Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)

        Int bam_to_fq_disk_size
        String bam_to_fq_mem_size    

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
    if (!is_bam_) {
        call cramToBam {
            input :
            cram_file = bam_or_cram_file,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker,
            disk_size = cram_to_bam_disk_size,
            mem_size = cram_to_bam_mem_size
        }
    }


    #################################################################################
    ####        Convert bam to fq1 & fq2                                            #
    #################################################################################
    call bamToFq { 
        input :
            bam_file = select_first([cramToBam.bam_file, bam_or_cram_file]),
            downsample_docker = downsample_docker,
            disk_size = bam_to_fq_disk_size,
            mem_size = bam_to_fq_mem_size
    }



    output {
        File? outputBam = cramToBam.bam_file
        File fastq_1 = bamToFq.fastq_file_1
        File fastq_2 = bamToFq.fastq_file_2
    }

}


task cramToBam {
    
    input {
        File cram_file
        File reference_fasta
        String downsample_docker
        Int? disk_size
        String? mem_size
    }

    File reference_index_file = reference_fasta + ".fai"
    
    String bam_file_name = basename(cram_file, ".cram") + ".bam"
    Int num_cpu = 4
    Float mem_size_gb = num_cpu * 4.0
  
    Float cram_inflate_ratio = 3.5
    Float disk_overhead = 20.0
    Float cram_size = size(cram_file, "GiB")
    Float bam_size = cram_inflate_ratio * cram_size
    Float ref_size = size(reference_fasta, "GiB")
    Float ref_index_size = size(reference_index_file, "GiB")
    Int vm_disk_size = ceil(bam_size + ref_size + ref_index_size + disk_overhead)

    output {
        File bam_file = bam_file_name
    }

    command {
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
    cpu: num_cpu
    memory: mem_size_gb
    disks: "local-disk " + vm_disk_size + " HDD"
    bootDiskSizeGb: 10
    preemptible: 3
    maxRetries: 1
    }
}


task bamToFq {
  
    input {
        File bam_file
        String downsample_docker
        Int disk_size
        String mem_size
    }

    String fastq_file_1_name = basename(bam_file, ".bam") + "_1.fastq"
    String fastq_file_2_name = basename(bam_file, ".bam") + "_2.fastq"

    output {
        File fastq_file_1 = fastq_file_1_name
        File fastq_file_2 = fastq_file_2_name
    }

    command <<<
        picard SamToFastq \
                I=~{bam_file} \
                FASTQ=~{fastq_file_1_name} \
                SECOND_END_FASTQ=~{fastq_file_2_name}
    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
    }
}



