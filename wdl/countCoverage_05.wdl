version 1.0

workflow countCoverage {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File downsample_sorted_cram
        File reference_fasta
        
        # Docker
        String downsample_docker

    }

    parameter_meta {
      downsample_sorted_cram: ".cram downsample and sorted file calculate WGS coverage (picard)."
      reference_fasta: ".fasta file with reference used to align bam or cram file"
    }

    meta {
        author: "Stephanie Hao"
        email: "shao@broadinstitute.org"
    }
    

    #################################################################################
    ####        Checks for downsamples sequences via picard                         #
    #################################################################################
    call countCoverage {
        input : 
            downsample_sorted_cram = downsample_sorted_cram,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker

    }

}

task countCoverage {
    
    input {
        File downsample_sorted_cram
        File reference_fasta
        String downsample_docker
    }

    String wgsCoverage_name = basename(downsample_sorted_cram, ".cram") + "_coverage.txt"

    output {
        File wgsCoverage = wgsCoverage_name
    }

    command {
        picard CollectWgsMetrics \
        I=${downsample_sorted_cram} \
        O=${wgsCoverage} \
        R=${reference_fasta} \
        COUNT_UNPAIRED=TRUE
    }

    runtime {
        docker: downsample_docker
    }
}