version 1.0

workflow bamToFq {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {

        File bam_file
        File reference_fasta
        # Docker
        String downsample_docker
        Int bam_to_fq_disk_size
        String bam_to_fq_mem_size

    }

    parameter_meta {
      bam_file: ".bam file to search for SVs. bams are preferable, crams will be converted to bams."
      reference_fasta: ".fasta file with reference used to align bam or cram file"
    }

    meta {
        author: "Stephanie Hao"
        email: "shao@broadinstitute.org"
    }
    

    #################################################################################
    ####        Convert bam to fq1 & fq2                                            #
    #################################################################################
    call bamToFq { 
        input :
            bam_file = bam_file,
            downsample_docker = downsample_docker,
            disk_size = bam_to_fq_disk_size,
            mem_size = bam_to_fq_mem_size
    }

    output {
        File fastq_1 = bamToFq.fastq_file_1
        File fastq_2 = bamToFq.fastq_file_2
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
