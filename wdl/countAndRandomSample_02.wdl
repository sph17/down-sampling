version 1.0

workflow countAndRandomSample {

    #################################################################################
    ####        Required basic arguments                                            #
    #################################################################################
    
    input {
        
        File reference_fasta
        File fastq_file_1
        File fastq_file_2
        Float start_depth
        Float final_depth
        Int? seed_override   
        # Docker
        String downsample_docker
        Int sampling_disk_size
        String sampling_mem_size
    }

    parameter_meta {
      fastq_file_: "paired .fastq file to downsample from."
      reference_fasta: ".fasta file with reference used to align bam or cram file"
      start_depth: "float/integer for initial depth of sequencing."
      final_depth: "float/integer for final desired depth."
    }

    meta {
        author: "Stephanie Hao"
        email: "shao@broadinstitute.org"
    }

    #################################################################################
    ####        Calculate reads and random samples appropriate # of reads           #
    #################################################################################
    call countAndRandomSample {
        input : 
            fastq_file_1 = fastq_file_1,
            fastq_file_2 = fastq_file_2,
            downsample_docker = downsample_docker,
            start_depth = start_depth,
            final_depth = final_depth,
            seed = select_first([seed_override, 20937]),
            disk_size = sampling_disk_size,
            mem_size = sampling_mem_size
    }

    output {
        File downsample_fastq_1 = countAndRandomSample.downsample_file_1
        File downsample_fastq_2 = countAndRandomSample.downsample_file_2
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
        File downsample_file_1 = fastq_downsample_1_name
        File downsample_file_2 = fastq_downsample_2_name
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

