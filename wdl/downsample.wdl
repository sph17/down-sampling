version 1.0

workflow downSampling {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File bam_or_cram_file
        File reference_fasta
        
        # Docker
        String downsample_docker

        Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
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
            downsample_docker = downsample_docker
        }
    }

    File bam_file = select_first([cramToBam.bam_file, bam_or_cram_file])


    #################################################################################
    ####        Convert bam to fq1 & fq2                                            #
    #################################################################################
    call bamToFq { 
        input :
            bam_file=bam_file,
            downsample_docker = downsample_docker
    }


    #################################################################################
    ####        Calculate reads and random samples appropriate # of reads           #
    #################################################################################
    call countAndRandomSample {
        input : 
            fastq_file_1=bamToFq.fastq_file_1,
            fastq_file_2=bamToFq.fastq_file_2,
            downsample_docker = downsample_docker
    }

    #################################################################################
    ####        bwa realigns the downsampled fastq files to Hg38 genome             #
    #################################################################################
    call realign {
        input :
            downsample_file_1 = countAndRandomSample.downsample_file_1,
            downsample_file_2 = countAndRandomSample.downsample_file_2,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker
    }

    #################################################################################
    ####        Sorts and indexes realigned sequences and convert to crams          #
    #################################################################################
    call sortIndex {
        input :
            cram_downsample_file = realign.cram_downsample_file,
            downsample_docker = downsample_docker

    }

    #################################################################################
    ####        Checks for downsamples sequences via picard                         #
    #################################################################################
    call countCoverage {
        input : 
            downsample_sorted_cram = sortIndex.cram_sorted_file,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker

    }

}


task cramToBam {
    
    input {
        File cram_file
        File reference_fasta
        String downsample_docker
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
    }
}

task bamToFq {
  
    input {
        File bam_file
        String downsample_docker
    }

    String fastq_file_1_name = basename(bam_file, ".bam") + "_1.fastq"
    String fastq_file_2_name = basename(bam_file, ".bam") + "_2.fastq"

    output {
        File fastq_file_1 = fastq_file_1_name
        File fastq_file_2 = fastq_file_2_name
    }

    command {
        picard SamToFastq \
                I=~{bam_file} \
                FASTQ=~{fastq_file_1_name} \
                SECOND_END_FASTQ=~{fastq_file_2_name}
    }

    runtime {
        docker: downsample_docker
    }
}

task countAndRandomSample {
    
    input {
        File fastq_file_1
        File fastq_file_2
        Float start_depth
        Float final_depth
        String downsample_docker
    }

    String fastq_downsample_1_name = basename(fastq_file_1, ".fastq") + "_downsample.fastq"
    String fastq_downsample_2_name = basename(fastq_file_2, ".fastq") + "_downsample.fastq"

    output {
        File downsample_file_1 = fastq_downsample_1_name
        File downsample_file_2 = fastq_downsample_2_name
    }

    command <<<

        count1=$(bash /opt/count_fastq.sh ${fastq_file_1})
        echo ${count1}
        count2=$(bash /opt/count_fastq.sh ${fastq_file_2})
        echo ${count2}

        initial=$(echo $start_depth | awk ' { printf "%0.2f\n", ($1 / 2); } ')
        final=$(echo $final_depth | awk ' { printf "%0.2f\n", ($1 / 2); } ')

        if [ $count1 -eq $count2 ]
          then

            freads=$(bash /opt/calcRD.sh $initial $final $count1)  #half of the total read coverage 30x

            echo "freads: " $freads

            seed=$RANDOM
            echo "seed: " $seed

            fastq-sample -n ${freads} --seed ${seed} -o ${downsample_file_1} ${fastq_file_1}

            fastq-sample -n ${freads} --seed ${seed} -o ${downsample_file_2} ${fastq_file_2}

        else
            echo "Error: counts don't match up!"
            echo $count1
            echo $count2
        fi
    >>>

    runtime {
        docker: downsample_docker
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
