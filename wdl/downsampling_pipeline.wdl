version 1.0

workflow downSampling {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        File bam_or_cram_file
        File reference_fasta
        Int cram_to_bam_disk_size
        String cram_to_bam_mem_size
        String downsample_docker

        Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)

        Int bam_to_fq_disk_size
        String bam_to_fq_mem_size

        Float start_depth
        Float final_depth
        Int? seed_override   
        Int sampling_disk_size
        String sampling_mem_size

        
        Int realign_disk_size
        String realign_mem_size

        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fai
        File ref_dict

        Int sort_disk_size
        String sort_mem_size

        Int coverage_disk_size
        String coverage_mem_size    

    }

    parameter_meta {
      bam_or_cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
      bam_file: ".bam file to search for SVs. bams are preferable, crams will be converted to bams."
      reference_fasta: ".fasta file with reference used to align bam or cram file"
      fastq_file_: "paired .fastq file to downsample from."
      start_depth: "float/integer for initial depth of sequencing."
      final_depth: "float/integer for final desired depth."
      cram_downsample_file: "randomly downsampled .cram files."
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


    #################################################################################
    ####        Calculate reads and random samples appropriate # of reads           #
    #################################################################################
    call countAndRandomSample {
        input : 
            fastq_file_1 = bamToFq.fastq_file_1,
            fastq_file_2 = bamToFq.fastq_file_2,
            downsample_docker = downsample_docker,
            start_depth = start_depth,
            final_depth = final_depth,
            seed = select_first([seed_override, 20937]),
            disk_size = sampling_disk_size,
            mem_size = sampling_mem_size
    }


    #################################################################################
    ####        bwa realigns the downsampled fastq files to Hg38 genome             #
    #################################################################################
    call realign {
        input :
            downsample_file_1 = countAndRandomSample.downsample_file_1,
            downsample_file_2 = countAndRandomSample.downsample_file_2,
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
            ref_dict = ref_dict

    }


    #################################################################################
    ####        Sorts and indexes realigned sequences and convert to crams          #
    #################################################################################
    call sortIndex {
        input :
            cram_downsample_file = realign.cram_downsample_file,
            downsample_docker = downsample_docker,
            disk_size = sort_disk_size,
            mem_size = sort_mem_size

    }


    #################################################################################
    ####        Checks for downsamples sequences via picard                         #
    #################################################################################
    call countCoverage {
        input : 
            downsample_sorted_cram = sortIndex.cram_sorted_file,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker,
            disk_size = coverage_disk_size,
            mem_size = coverage_mem_size,
            reference_dict = ref_dict,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            ref_fai = ref_fai 

    }

    output {
        File? outputBam = cramToBam.bam_file
        File fastq_1 = bamToFq.fastq_file_1
        File fastq_2 = bamToFq.fastq_file_2
        File downsample_fastq_1 = countAndRandomSample.downsample_file_1
        File downsample_fastq_2 = countAndRandomSample.downsample_file_2
        File downsample_cram = realign.cram_downsample_file
        File sorted_cram = sortIndex.cram_sorted_file
        File crai_file = sortIndex.crai_file
        File wgsCoverage_metrics = countCoverage.wgsCoverage
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

task countAndRandomSample {
    
    input {
        File fastq_file_1
        File fastq_file_2
        Float start_depth
        Float final_depth
        String downsample_docker
        Int seed
        Int disk_size
        String mem_size
    }

    String fastq_downsample_1_name = basename(fastq_file_1, ".fastq") + "_downsample"
    String fastq_downsample_2_name = basename(fastq_file_2, ".fastq") + "_downsample"


    output {
        File downsample_file_1 = fastq_downsample_1_name + ".fastq"
        File downsample_file_2 = fastq_downsample_2_name + ".fastq"
    }

    command <<<

        count1=$(bash /opt/count_fastq.sh ~{fastq_file_1})
        echo ${count1}
        count2=$(bash /opt/count_fastq.sh ~{fastq_file_2})
        echo ${count2}

        initial=$(echo ~{start_depth} | awk ' { printf "%0.2f\n", ($1 / 2); } ')
        final=$(echo ~{final_depth} | awk ' { printf "%0.2f\n", ($1 / 2); } ')
        echo ${initial}
        echo ${final}

        if [ $count1 -eq $count2 ]
          then

            freads=$(bash /opt/calcRD.sh ${initial} ${final} ${count1})  #half of the total read coverage 30x

            echo "freads: ${freads}"


            fastq-sample -n ${freads} --seed ~{seed} -o ~{fastq_downsample_1_name} ~{fastq_file_1}

            fastq-sample -n ${freads} --seed ~{seed} -o ~{fastq_downsample_2_name} ~{fastq_file_2}

        else
            echo "Error: counts don't match up!"
            echo ${count1}
            echo ${count2}
        fi
    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
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
    }
    
    String bam_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample.bam"
    String cram_downsample_name = basename(downsample_file_1, "_1_downsample.fastq") + "_downsample.cram"

    
    output {
        File cram_downsample_file = cram_downsample_name
    }

    command <<<
        NT=$(nproc)

        bwa mem -M \
                -t ${NT}\
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
        cpu: "5"
        disks: "local-disk " + disk_size + " HDD"
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
        File crai_file = downsample_file_sorted_name + ".crai"
    }

    command <<<

        samtools sort -o ~{downsample_file_sorted_name} ~{cram_downsample_file}

        samtools index ~{downsample_file_sorted_name}

    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
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



