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
    }

    parameter_meta {
      downsample_file_: "downsampled .fastq files to be realgined"
      reference_fasta: ".fasta file with reference used to align bam or cram file"
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
            downsample_docker = downsample_docker
    }
}

task realign {
    
    input {
        File downsample_file_1 
        File downsample_file_2 
        File reference_fasta
        String downsample_docker
    }
    
    String bam_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample.bam"
    String cram_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample.cram"

    
    output {
        File cram_downsample_file = cram_downsample_name
    }

    command <<<
        NT=$(nproc)

        bwa mem -C \
                -t ${NT} \
                ${reference_fasta} \
                ${downsample_file_1} \
                ${downsample_file_2} \
                | samtools view \
                -bS - > ${bam_downsample_name}

        samtools view \
                -C \
                -T ${reference_fasta} \
                ${bam_downsample_name} \
                > ${cram_downsample_name}
    >>>

    runtime {
        docker: downsample_docker
    }
}