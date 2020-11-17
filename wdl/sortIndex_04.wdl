version 1.0

workflow sortIndex {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File cram_downsample_file
        File reference_fasta
        Int sort_disk_size
        String sort_mem_size
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
            downsample_docker = downsample_docker,
            disk_size = sort_disk_size,
            mem_size = sort_mem_size

    }

    output {
        File sorted_cram = sortIndex.cram_sorted_file
        File crai_file = sortIndex.crai_file
    }

}

task sortIndex {
    
    input {

        File cram_downsample_file
        String downsample_docker
        Int disk_size
        String mem_size
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
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
    }
}
