version 1.0

workflow downSampling_02 {

    #################################################################################
    ####        Required basic arguments for downsampling part 2                    #
    #################################################################################
    
    input {
        File reference_fasta
        String downsample_docker
        File fastq_file_1
        File fastq_file_2
        Float start_depth
        Float final_depth
        Int? seed_override   
        Int sampling_disk_size
        String sampling_mem_size
        File original_cram_file
        
        Int realign_disk_size
        String realign_mem_size

        Int addRG_disk_size
        String addRG_mem_size

        Int markDup_disk_size
        String markDup_mem_size

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

        String sample_ID
        File intervals_exons

        File? gatk4_jar_override
        String gatk_docker  

    }

    parameter_meta {
      original_cram_file: "The original cram file needed for re-adding RG."
      reference_fasta: ".fasta file with reference used to align bam or cram file"
      fastq_file_: "paired .fastq file to downsample from."
      start_depth: "float/integer for initial depth of sequencing."
      final_depth: "float/integer for final desired depth."
      seed: "optional seed input for random sampling."
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


    #################################################################################
    ####        bwa realigns the downsampled fastq files to Hg38 genome,            #
    ####        adds read groups and mark duplicates and converts to cram           #
    #################################################################################
    call realign {
        input :
            downsample_file_1 = countAndRandomSample.downsample_file_1,
            downsample_file_2 = countAndRandomSample.downsample_file_2,
            reference_fasta = reference_fasta,
            downsample_docker = downsample_docker,
            disk_size = realign_disk_size,
            mem_size = realign_mem_size,
            final_depth = final_depth,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            ref_fai = ref_fai,
            ref_dict = ref_dict
    }

    call addRG {
        input :
            bam_downsample_file = realign.bam_downsample_file,
            downsample_docker = downsample_docker,
            disk_size = addRG_disk_size,
            mem_size = addRG_mem_size,
            original_cram_file = original_cram_file

    }

    call markDuplicates {
        input :
            bam_sorted_rg_file = addRG.bam_sorted_rg_file,
            downsample_docker = downsample_docker,
            disk_size = markDup_disk_size,
            mem_size = markDup_mem_size,
            reference_fasta = reference_fasta

    }


    #################################################################################
    ####        Sorts and indexes realigned sequences and convert to crams          #
    #################################################################################
    call sortIndex {
        input :
            cram_downsample_file = markDuplicates.cram_downsample_file,
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

    #################################################################################
    ####        Collect Coverage with gatk                                          #
    #################################################################################
    call CollectCountsCram {
        input :
            intervals_exons = intervals_exons,
            cram = sortIndex.cram_sorted_file,
            crai = sortIndex.crai_file,
            sample_ID = sample_ID,
            hg38_reference = reference_fasta,
            hg38_reference_fai = ref_fai,
            hg38_reference_dict = ref_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker

    }


    output {
        File downsample_fastq_1 = countAndRandomSample.downsample_file_1
        File downsample_fastq_2 = countAndRandomSample.downsample_file_2
        File downsample_cram = markDuplicates.cram_downsample_file
        File markdup_metrics = markDuplicates.markdup_metrics_file
        File sorted_cram = sortIndex.cram_sorted_file
        File crai_file = sortIndex.crai_file
        File wgsCoverage_metrics = countCoverage.wgsCoverage
        File read_counts = CollectCountsCram.counts_reads
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

    String fastq_downsample_1_name = basename(fastq_file_1, ".fastq") + "_downsample_~{final_depth}x"
    String fastq_downsample_2_name = basename(fastq_file_2, ".fastq") + "_downsample_~{final_depth}x"


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
        Float final_depth
    }
    
    String bam_downsample_name = basename(downsample_file_1, "_1_downsample_~{final_depth}x.fastq") + "_downsample.bam"

    
    output {
        File bam_downsample_file = bam_downsample_name
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
        
    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "5"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task addRG {
    
    input {
        File bam_downsample_file
        String downsample_docker
        Int disk_size
        String mem_size
        File original_cram_file
    }
    
    String bam_readgroup_name = basename(bam_downsample_file, "_downsample.bam") + "_downsample_rg.bam"
    String bam_sorted_name = basename(bam_downsample_file, "_downsample.bam") + "_downsample_rg_sorted.bam"

    
    output {
        File bam_sorted_rg_file = bam_sorted_name
    }

    command <<<
       
        #re-add read group
        SM=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f7 | sed 's/SM://g')

        PU=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f9 | sed 's/PU://g')

        ID=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f2 | sed 's/ID://g')

        LB=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f5 | sed 's/LB://g')

        PL=$(samtools view -H ~{original_cram_file} | grep '^@RG' | head -1 | cut -f3 | sed 's/PL://g')

        picard AddOrReplaceReadGroups \
            I=~{bam_downsample_file} \
            O=~{bam_readgroup_name} \
            RGID=${ID} \
            RGLB=${LB} \
            RGPL=${PL} \
            RGPU=${PU} \
            RGSM=${SM}

        #sort bam for mark dup
        samtools sort -o ~{bam_sorted_name} ~{bam_readgroup_name}

    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        maxRetries: 1
    }
}

task markDuplicates {
    
    input {
        File bam_sorted_rg_file
        String downsample_docker
        Int disk_size
        String mem_size
        File reference_fasta
    }
    
    String bam_markdup_name = basename(bam_sorted_rg_file, "_downsample_rg_sorted.bam") + "_downsample_rg_sorted_markdup.bam"
    String markdup_metrics_name = basename(bam_sorted_rg_file, "_downsample_rg_sorted.bam") + "_downsample_markdup_metrics.txt"
    String cram_downsample_name = basename(bam_sorted_rg_file, "_downsample_rg_sorted.bam") + "_downsample_markdup.cram"
    
    output {
        File cram_downsample_file = cram_downsample_name
        File markdup_metrics_file = markdup_metrics_name
    }

    command <<<

        #marking duplicates
        picard MarkDuplicates \
            I=~{bam_sorted_rg_file} \
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
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        maxRetries: 1
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

        echo "File: sorted"

        samtools index ~{downsample_file_sorted_name}

        echo "File: indexed"        

    >>>

    runtime {
        docker: downsample_docker
        memory: mem_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        maxRetries: 1
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
        preemptible: 3
        maxRetries: 1
    }
}

task CollectCountsCram {
    input {
        File intervals_exons
        File cram
        File crai
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        File? gatk4_jar_override 
        String sample_ID
        String gatk_docker
        Int? disk_space_gb
    }

    # Runtime parameters

    Boolean use_ssd = false
    Int num_cpu = 1
    Int machine_mem_gb = 3
    Int command_mem_gb = machine_mem_gb
    String base_filename = basename(cram, ".cram")
    String counts_reads_filename = "${sample_ID}.counts.tsv"
    
    command <<<
        set -euo pipefail
        export GATK_LOCAL_JAR="/root/gatk.jar"

        gatk --java-options "-Xmx~{command_mem_gb}G" CollectReadCounts \
            -I ~{cram} \
            --read-index ~{crai} \
            -L ~{intervals_exons} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --reference ~{hg38_reference} \
            --format TSV \
            -O ~{counts_reads_filename}

    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_gb + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(cram, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: num_cpu
        preemptible: 3
        maxRetries: 1
    }

    output {
        File counts_reads = counts_reads_filename
    }
}



