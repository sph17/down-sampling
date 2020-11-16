version 1.0

workflow cramToBam {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File bam_or_cram_file
        File reference_fasta
        
        # Docker
        String downsample_docker

        Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
    	Int cram_to_bam_disk_size
		String cram_to_bam_mem_size
    }

    parameter_meta {
      cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
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

    output {
    	File? outputBam = cramToBam.bam_file
	}

}


task cramToBam {
    
    input {
        File cram_file
        File reference_fasta
        String downsample_docker
        Int disk_size
		String mem_size

    }

    String bam_file_name = basename(cram_file, ".cram") + ".bam"

    output {
        File bam_file = bam_file_name
    }

    command {
        samtools view \
                -b \
                -o ~{bam_file_name} \
                ~{cram_file}
    }

    runtime {
    docker: downsample_docker
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
	}
}

