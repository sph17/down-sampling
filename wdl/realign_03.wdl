version 1.0

workflow realign {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File reference_fasta
        File downsample_file_1
        File downsample_file_2
        File original_cram_file
        # Docker
        String downsample_docker
        Int realign_disk_size
        String realign_mem_size

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
    ####        bwa realigns the downsampled fastq files to Hg38 genome,            #
    ####        adds read groups and mark duplicates and converts to cram           #
    #################################################################################
    call realign {
        input :
            downsample_file_1 = downsample_file_1,
            downsample_file_2 = downsample_file_2,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker,
            disk_size = realign_disk_size,
            mem_size = realign_mem_size,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            original_cram_file = original_cram_file

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
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fai
        File ref_dict
        File original_cram_file
    }
    
    String bam_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample.bam"
    String bam_readgroup_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample_rg.bam"
    String bam_sorted_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample_rg_sorted.bam"
    String bam_markdup_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample_rg_sorted_markdup.bam"
    String markdup_metrics_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample_markdup_metrics.txt"
    String cram_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample_markdup.cram"
    
    output {
        File cram_downsample_file = cram_downsample_name
        File markdup_metrics_file = markdup_metrics_name
    }

    command <<<
        #realign
        NT=$(nproc)

        bwa mem -M \
                -t ${NT}\
                ~{reference_fasta} \
                ~{downsample_file_1} \
                ~{downsample_file_2} \
                | samtools view \
                -bS - > ~{bam_downsample_name}
        

        #re-add read group
        SM=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f7 | sed 's/SM://g')

        PU=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f9 | sed 's/PU://g')

        ID=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f2 | sed 's/ID://g')

        LB=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f5 | sed 's/LB://g')

        PL=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f3 | sed 's/PL://g')

        picard AddOrReplaceReadGroups \
            I=~{bam_downsample_name} \
            O=~{bam_readgroup_name} \
            RGID=${ID} \
            RGLB=${LB} \
            RGPL=${PL} \
            RGPU=${PU} \
            RGSM=${SM}

        #sort bam for mark dup
        samtools sort -o ~{bam_sorted_name} ~{bam_readgroup_name}

        #marking duplicates
        picard MarkDuplicates \
            I=~{bam_sorted_name} \
            O=~{bam_markdup_name} \
            M=~{markdup_metrics_name}

        #bam to cram
        samtools view \
                -C \
                -T ~{reference_fasta} \
                ~{bam_markdup_name} \
                > ~{cram_downsample_name}
    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "5"
        disks: "local-disk " + disk_size + " HDD"
    }
}