version 1.0

workflow sortIndex {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File cram_downsample_file
        File reference_fasta
        
        # Docker
        String downsample_docker
    }

    parameter_meta {
      cram_downsample_file: "randomly downsampled .cram files."
      reference_fasta: ".fasta file with reference used to align bam or cram file"
    }

    meta {
        author: "Stephanie Hao"
        email: "shao@broadinstitute.org"
    }
    

    #################################################################################
    ####        Sorts and indexes realigned sequences and convert to crams          #
    #################################################################################
    call sortIndex {
        input :
            cram_downsample_file = cram_downsample_file,
            downsample_docker = downsample_docker

    }

}

task sortIndex {
    
    input {

        File cram_downsample_file
        String downsample_docker
    }

    String downsample_file_sorted_name = basename(cram_downsample_file, ".cram") + "_sorted.cram"
    
    output {
        File cram_sorted_file = downsample_file_sorted_name
        File crai_file = cram_sorted_file + ".crai"
    }

    command <<<

        samtools sort -o ${downsample_file_sorted_name} ${cram_downsample_file}

        samtools index ${downsample_file_sorted_name}

    >>>

    runtime {
        docker: downsample_docker
    }
}
