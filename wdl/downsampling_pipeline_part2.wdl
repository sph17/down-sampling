version 1.0

import "Structs.wdl"

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
    File original_cram_or_bam_file_read_groups

    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai
    File ref_dict

    String sample_ID
    File intervals_genome

    File? gatk4_jar_override
    String gatk_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_random_sample
    RuntimeAttr? runtime_attr_realign
    RuntimeAttr? runtime_attr_add_read_group
    RuntimeAttr? runtime_attr_mark_duplicates
    RuntimeAttr? runtime_attr_sort_index
    RuntimeAttr? runtime_attr_count_coverage
    RuntimeAttr? runtime_attr_collect_counts
    RuntimeAttr? runtime_attr_depth_of_coverage

    Boolean run_count_coverage = true
    Boolean run_depth_of_coverage = false

  }

  parameter_meta {
    original_cram_or_bam_file_read_groups: "The original cram file's read group extracted, it is needed to re-add RG."
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
      runtime_attr_override = runtime_attr_random_sample
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
      runtime_attr_override = runtime_attr_realign,
      final_depth = final_depth,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      ref_fai = ref_fai,
      ref_dict = ref_dict
  }

  #################################################################################
  ####        Re-adds Read Groups back to realigned bam files                     #
  #################################################################################
  call addReadGroupAndSort {
    input :
      bam_downsample_file = realign.bam_downsample_file,
      downsample_docker = downsample_docker,
      reference_fasta = reference_fasta,
      original_cram_or_bam_file_read_groups = original_cram_or_bam_file_read_groups,
      runtime_attr_override = runtime_attr_add_read_group
  }

  #################################################################################
  ####        Marks duplicates to closely mimick format of sequencers             #
  ####        and converts to cram                                                #
  #################################################################################
  call markDuplicatesAndToCram {
    input :
      bam_sorted_rg_file = addReadGroupAndSort.bam_sorted_rg_file,
      downsample_docker = downsample_docker,
      runtime_attr_override = runtime_attr_mark_duplicates,
      reference_fasta = reference_fasta
  }


  #################################################################################
  ####        Sorts and indexes realigned sequences                               #
  #################################################################################
  call sortIndex {
    input :
      cram_downsample_file = markDuplicatesAndToCram.cram_downsample_file,
      downsample_docker = downsample_docker,
      runtime_attr_override = runtime_attr_sort_index
  }


  #################################################################################
  ####        Checks for downsamples sequences via picard                         #
  #################################################################################
  if (run_count_coverage){
    call countCoverage {
    input : 
      downsample_sorted_cram = sortIndex.cram_sorted_file,
      reference_fasta = reference_fasta,
      downsample_docker = downsample_docker,
      runtime_attr_override = runtime_attr_count_coverage,
      reference_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      ref_fai = ref_fai 
    }
  }
  

  #################################################################################
  ####        Collect Read Counts with gatk                                       #
  #################################################################################
  call collectCountsCram {
    input :
      intervals_genome = intervals_genome,
      cram = sortIndex.cram_sorted_file,
      crai = sortIndex.crai_file,
      sample_ID = sample_ID,
      hg38_reference = reference_fasta,
      hg38_reference_fai = ref_fai,
      hg38_reference_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_collect_counts
  }

  #################################################################################
  ####        Calculates per locus depth of coverage with gatk                    #
  #################################################################################
  if (run_depth_of_coverage){
    call calculateDepthOfCoverage {
    input : 
      downsample_sorted_cram = sortIndex.cram_sorted_file,
      crai = sortIndex.crai_file,
      reference_fasta = reference_fasta,
      gatk_docker = gatk_docker,
      sample_ID = sample_ID,
      intervals_genome = intervals_genome,
      runtime_attr_override = runtime_attr_depth_of_coverage,
      reference_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      ref_fai = ref_fai
    }
  }


  output {
    File markdup_metrics = markDuplicatesAndToCram.markdup_metrics_file
    File sorted_cram = sortIndex.cram_sorted_file
    File crai_file = sortIndex.crai_file
    File? wgs_coverage_metrics = countCoverage.wgs_coverage_file
    File read_counts = collectCountsCram.counts_reads
    File? depth_of_coverage = calculateDepthOfCoverage.depth_of_coverage_file

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
    RuntimeAttr? runtime_attr_override
  }

  Int num_cpu = 1
  Int mem_size_gb = 6
  Int vm_disk_size = 400

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }  

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String fastq_downsample_1_name = basename(fastq_file_1, "_1.fastq") + "_downsample_~{final_depth}x"
  String fastq_downsample_2_name = basename(fastq_file_2, "_2.fastq") + "_downsample_~{final_depth}x"

  output {
    File downsample_file_1 = fastq_downsample_1_name + ".1.fastq"
    File downsample_file_2 = fastq_downsample_2_name + ".2.fastq"
  }

  command <<<
    set -euo pipefail

    #counts the number of reads in each fastq file and prints for qc
    count1=$(bash /opt/count_fastq.sh ~{fastq_file_1})
    echo ${count1} 

    count2=$(bash /opt/count_fastq.sh ~{fastq_file_2})
    echo ${count2}

    #calculates half of the initial and final depth for the number of reads to sample in each fq file
    initial=$(echo ~{start_depth} | awk ' { printf "%0.2f\n", ($1 / 2); } ')
    final=$(echo ~{final_depth} | awk ' { printf "%0.2f\n", ($1 / 2); } ')
    echo ${initial}
    echo ${final}

    if [ $count1 -eq $count2 ] #checks that the paired-end fastq files have the same number of reads
      then
      #half of the total read coverage (i.e. 15x for 30x)
        freads=$(bash /opt/calcRD.sh ${initial} ${final} ${count1})  

        echo "freads: ${freads}"

        #uses now the calculated freads (1/2 of final read coverage wanted) in each fastq file
        fastq-sample -n ${freads} --seed ~{seed} -o ~{fastq_downsample_1_name} ~{fastq_file_1} ~{fastq_file_2}
        #output name will be: sample_id.final_downsample_~{final_depth}x.[1/2].fastq

      else
        echo "Error: counts don't match up!"
        echo ${count1}
        echo ${count2}
    fi
  >>>

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task realign {
    
  input {
    File downsample_file_1 
    File downsample_file_2 
    File reference_fasta
    String downsample_docker
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai
    File ref_dict
    Float final_depth
    RuntimeAttr? runtime_attr_override
  }
  
  Int num_cpu = 16
  Int mem_size_gb = 21
  Int vm_disk_size = 100

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }  

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String bam_downsample_name = basename(downsample_file_1, "_downsample_~{final_depth}x.1.fastq") + "_downsample.bam"
  
  output {
    File bam_downsample_file = bam_downsample_name
  }

  command <<<
    set -euo pipefail

    #realigns downsampled paired fq files to sam and then bam with bwa
    NT=$(nproc)

    bwa mem -K 100000000 -Y -M \
      -t ${NT}\
      ~{reference_fasta} \
      ~{downsample_file_1} \
      ~{downsample_file_2} \
      | samtools view \
      -bS - > ~{bam_downsample_name}
        
  >>>

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task addReadGroupAndSort {
    
  input {
    File bam_downsample_file
    String downsample_docker
    File reference_fasta
    File original_cram_or_bam_file_read_groups
    RuntimeAttr? runtime_attr_override
  }
  
  Int num_cpu = 1
  Int mem_size_gb = 3
  Int vm_disk_size = 70

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }  

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String bam_readgroup_name = basename(bam_downsample_file, "_downsample.bam") + "_downsample_rg.bam"
  String bam_sorted_name = basename(bam_downsample_file, "_downsample.bam") + "_downsample_rg_sorted.bam"

    
  output {
    File bam_sorted_rg_file = bam_sorted_name
  }

  command <<<
    set -euo pipefail
    
    #Re-add read group, necessary for collectCountCrams. This is not best practice. This only adds one RG line, when there are multiples
    SM=$(grep '^@RG' ~{original_cram_or_bam_file_read_groups} | head -1 | cut -f7 | sed 's/SM://g')

    PU=$(grep '^@RG' ~{original_cram_or_bam_file_read_groups} | head -1 | cut -f9 | sed 's/PU://g')

    ID=$(grep '^@RG' ~{original_cram_or_bam_file_read_groups} | head -1 | cut -f2 | sed 's/ID://g')

    LB=$(grep '^@RG' ~{original_cram_or_bam_file_read_groups} | head -1 | cut -f5 | sed 's/LB://g')

    PL=$(grep '^@RG' ~{original_cram_or_bam_file_read_groups} | head -1 | cut -f3 | sed 's/PL://g')

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
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task markDuplicatesAndToCram {
  
  input {
    File bam_sorted_rg_file
    String downsample_docker
    File reference_fasta
    RuntimeAttr? runtime_attr_override
  }
  
  Int num_cpu = 1
  Int mem_size_gb = 3
  Int vm_disk_size = 50

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String bam_markdup_name = basename(bam_sorted_rg_file, "_downsample_rg_sorted.bam") + "_downsample_rg_sorted_markdup.bam"
  String markdup_metrics_name = basename(bam_sorted_rg_file, "_downsample_rg_sorted.bam") + "_downsample_markdup_metrics.txt"
  String cram_downsample_name = basename(bam_sorted_rg_file, "_downsample_rg_sorted.bam") + "_downsample_markdup.cram"
    
  output {
    File cram_downsample_file = cram_downsample_name
    File markdup_metrics_file = markdup_metrics_name
  }

  command <<<
    set -euo pipefail
    #Identifies duplicate reads.
    picard MarkDuplicates \
      I=~{bam_sorted_rg_file} \
      O=~{bam_markdup_name} \
      M=~{markdup_metrics_name} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER=queryname

    #bam to cram
    samtools view \
      -C \
      -T ~{reference_fasta} \
      ~{bam_markdup_name} \
      > ~{cram_downsample_name}
  >>>

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task sortIndex {
    
  input {
    File cram_downsample_file
    String downsample_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_cpu = 1
  Int mem_size_gb = 4
  Int vm_disk_size = 50

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 20,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String downsample_file_sorted_name = basename(cram_downsample_file, ".cram") + "_sorted.cram"
    
  output {
    File cram_sorted_file = downsample_file_sorted_name
    File crai_file = downsample_file_sorted_name + ".crai"
  }

  command <<<
      
    set -euo pipefail
      
    #sorts and indexes the downsampled, markduped cram file

    samtools sort -o ~{downsample_file_sorted_name} ~{cram_downsample_file}

    echo "File: sorted"

    samtools index ~{downsample_file_sorted_name}

    echo "File: indexed"        

  >>>

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task countCoverage {
    
  input {
    File downsample_sorted_cram
    File reference_fasta
    String downsample_docker
    File reference_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai
    RuntimeAttr? runtime_attr_override
  }

  Int num_cpu = 1
  Int mem_size_gb = 4
  Int vm_disk_size = 50

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String wgs_coverage_name = basename(downsample_sorted_cram, ".cram") + "_coverage.txt"
  String ref_dictionary_name = "Homo_sapiens_assembly38.dict"

  output {
    File wgs_coverage_file = wgs_coverage_name
  }

  command {
    set -euo pipefail

    java -Xmx4G -jar /opt/conda/share/picard-2.23.8-0/picard.jar CollectWgsMetrics \
      I=~{downsample_sorted_cram} \
      O=~{wgs_coverage_name} \
      R=~{reference_fasta} \
      COUNT_UNPAIRED=TRUE
    }

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task collectCountsCram {
  input {
    File intervals_genome
    File cram
    File crai
    File hg38_reference
    File hg38_reference_fai
    File hg38_reference_dict
    String sample_ID
    String gatk_docker
    File? gatk4_jar_override 
    RuntimeAttr? runtime_attr_override
  }

  # Runtime parameters adapted from gatk-sv "CollectCoverage.wdl"
  Int num_cpu = 1
  Int mem_size_gb = 6
  Int vm_disk_size = 50

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_overhead_gb = 2.0
  Int command_mem_mb = floor((mem_size_gb - mem_overhead_gb) * 1024)


  String base_filename = basename(cram, ".cram")
  String counts_reads_filename = "${sample_ID}.counts.tsv"
    
  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR="/root/gatk.jar"
    #Collects read counts and output as a counts.tsv file for downstream analyses
    gatk --java-options "-Xmx~{command_mem_mb}m" CollectReadCounts \
      -I ~{cram} \
      --read-index ~{crai} \
      -L ~{intervals_genome} \
      --interval-merging-rule OVERLAPPING_ONLY \
      --reference ~{hg38_reference} \
      --format TSV \
      -O ~{counts_reads_filename}

  >>>

  runtime {
    docker: gatk_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File counts_reads = counts_reads_filename
  }
}

task calculateDepthOfCoverage {
  input {
    File downsample_sorted_cram
    File crai
    File reference_fasta
    String gatk_docker
    File intervals_genome
    File reference_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai
    String sample_ID
    RuntimeAttr? runtime_attr_override
  }

  Int num_cpu = 16
  Int mem_size_gb = 30
  Int vm_disk_size = 400

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String depth_coverage_name = "${sample_ID}_depth_of_coverage.txt"


  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR="/root/gatk.jar"
    #run on generated downsampled BAMs to produce per locus coverage
    gatk DepthOfCoverage \
      -I ~{downsample_sorted_cram} \
      -L ~{intervals_genome} \
      -O ~{depth_coverage_name} \
      -R ~{reference_fasta}


  >>>

  runtime {
    docker: gatk_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File depth_of_coverage_file = depth_coverage_name
  }

}

