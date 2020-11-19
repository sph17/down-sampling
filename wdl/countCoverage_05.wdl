version 1.0

workflow countCoverage {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File downsample_sorted_cram
        File reference_fasta
        Int coverage_disk_size
        String coverage_mem_size        
        # Docker
        String downsample_docker
        File reference_dict
        File ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fai


    }

    parameter_meta {
      downsample_sorted_cram: {help: ".cram downsample and sorted file calculate WGS coverage (picard)."}
      reference_fasta: { localization_optional: true, help: ".fasta file with reference used to align bam or cram file"}
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
            downsample_docker = downsample_docker,
            disk_size = coverage_disk_size,
            mem_size = coverage_mem_size,
            reference_dict = reference_dict,
            ref_alt = ref_alt,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sac = ref_sac,
            ref_fai = ref_fai 

    }

    output {
        File wgsCoverage_metrics = countCoverage.wgsCoverage
    }

}

task countCoverage {
    
    input {
        File downsample_sorted_cram
        File reference_fasta
        String downsample_docker
        Int disk_size
        String mem_size
        File reference_dict
        File ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fai
    }

    String wgsCoverage_name = basename(downsample_sorted_cram, ".cram") + "_coverage.txt"
    String ref_dictionary_name = "Homo_sapiens_assembly38.dict"

    output {
        File wgsCoverage = wgsCoverage_name
    }

    command {
        java -Xmx32G -jar /opt/conda/share/picard-2.23.8-0/picard.jar CollectWgsMetrics \
        I=~{downsample_sorted_cram} \
        O=~{wgsCoverage_name} \
        R=~{reference_fasta} \
        COUNT_UNPAIRED=TRUE
    }

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
    }
}
