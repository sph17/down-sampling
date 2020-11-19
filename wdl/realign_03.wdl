version 1.0

workflow realign {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File reference_fasta
        File downsample_file_1
        File downsample_file_2
        
        # Docker
        String downsample_docker
        Int realign_disk_size
        String realign_mem_size

        File ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fai
        File ref_dict
    }

    parameter_meta {
      downsample_file_: { help: "downsampled .fastq files to be realgined" }
      reference_fasta: { localization_optional: true, help: ".fasta file with reference used to align bam or cram file"}
    }

    meta {
        author: "Stephanie Hao"
        email: "shao@broadinstitute.org"
    }


    #################################################################################
    ####        bwa realigns the downsampled fastq files to Hg38 genome             #
    #################################################################################
    call realign {
        input :
            downsample_file_1 = downsample_file_1,
            downsample_file_2 = downsample_file_2,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker,
            disk_size = realign_disk_size,
            mem_size = realign_mem_size,
            ref_alt = ref_alt,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sac = ref_sac,
            ref_fai = ref_fai,
            ref_dict = ref_dict

    }

    output {
        File downsample_cram = realign.cram_downsample_file
    }
}

task realign {
    
    input {
        File downsample_file_1 
        File downsample_file_2 
        File reference_fasta
        String downsample_docker
        Int disk_size
        String mem_size
        File ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fai
        File ref_dict
    }
    
    String bam_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample.bam"
    String cram_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample.cram"

    
    output {
        File cram_downsample_file = cram_downsample_name
    }

    command <<<
        bwa mem -M \
                ~{reference_fasta} \
                ~{downsample_file_1} \
                ~{downsample_file_2} \
                | samtools view \
                -bS - > ~{bam_downsample_name}

        samtools view \
                -C \
                -T ~{reference_fasta} \
                ~{bam_downsample_name} \
                > ~{cram_downsample_name}
    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
    }
}