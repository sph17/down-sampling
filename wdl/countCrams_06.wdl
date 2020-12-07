version 1.0

workflow CC_Crams{
    input{

        Array[File] crams
        Array[File] crais
        Array[String] sample_IDs
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
    
        File intervals_exons

        File? gatk4_jar_override
        String gatk_docker
    }
    

    scatter (scatter_index in range(length(crams))){
        call CollectCountsCram {
            input :
                intervals_exons = intervals_exons,
                cram = crams[scatter_index],
                crai = crais[scatter_index],
				sample_ID = sample_IDs[scatter_index],
                hg38_reference = hg38_reference,
                hg38_reference_fai = hg38_reference_fai,
                hg38_reference_dict = hg38_reference_dict,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker

        }
    }

    output {
        Array [File] counts_exons = CollectCountsCram.counts_exons
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
    Int machine_mem_mb = 600
    Int command_mem_mb = machine_mem_mb - 100
    String base_filename = basename(cram, ".cram")
    String counts_exons_filename = "${sample_ID}.exons.counts.tsv"
	
    command <<<
        set -euo pipefail
        export GATK_LOCAL_JAR="/root/gatk.jar"

        gatk --java-options "-Xmx~{command_mem_mb}m" CollectReadCounts \
            -I ~{cram} \
            --read-index ~{crai} \
            -L ~{intervals_exons} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --reference ~{hg38_reference} \
            --format TSV \
            -O ~{counts_exons_filename}

    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(cram, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: num_cpu
    }

    output {
        File counts_exons = counts_exons_filename
    }
}


