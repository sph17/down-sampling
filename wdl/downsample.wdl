version 1.0

import "Structs.wdl"
import "downsampling_pipeline_part1.wdl" as ds1
import "downsampling_pipeline_part2.wdl" as ds2

workflow downSampling {

  #################################################################################
  ####        Required basic arguments for downsampling pipeline                  #
  #################################################################################
  
  input {
    File reference_fasta
    String downsample_docker
    File bam_or_cram_file

    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai
    File ref_dict

    String sample_ID
    File intervals_exons

    #default start depth is 30x
    Float start_depth = 30
    Int? seed_override  

    File? gatk4_jar_override
    String gatk_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_bam_to_fastq 
    RuntimeAttr? runtime_attr_random_sample
    RuntimeAttr? runtime_attr_realign
    RuntimeAttr? runtime_attr_add_read_group
    RuntimeAttr? runtime_attr_mark_duplicates
    RuntimeAttr? runtime_attr_sort_index
    RuntimeAttr? runtime_attr_count_coverage
    RuntimeAttr? runtime_attr_collect_counts   

    #Execution defaults and overrides
    Boolean run_downsample_2x = true
    Boolean run_downsample_4x = true
    Boolean run_downsample_6x = true
    Boolean run_downsample_8x = true
    Boolean run_downsample_custom = false

    Float? final_depth_custom

  }

  meta {
    author: "Stephanie Hao"
    email: "shao@broadinstitute.org"
  }
    
  #################################################################################
  ####        Calls part 1 to convert cram/bam file to paired fastq               #
  #################################################################################

  call ds1.downSampling_01 {
    input :
      bam_or_cram_file = bam_or_cram_file,
      reference_fasta = reference_fasta,
      downsample_docker = downsample_docker,
      runtime_attr_cram_to_bam = runtime_attr_cram_to_bam,
      runtime_attr_bam_to_fastq = runtime_attr_bam_to_fastq 
  }

  Float final_depth_2x = 2
  Float final_depth_4x = 4
  Float final_depth_6x = 6
  Float final_depth_8x = 8

  #################################################################################
  ####        Downsamples samples to 2x, 4x, 6x, and 8x read depths               #
  #################################################################################
  if (run_downsample_2x) {
    call ds2.downSampling_02 as downSampling_02_2x {
      input :
        fastq_file_1 = downSampling_01.fastq_1,
        fastq_file_2 = downSampling_01.fastq_2,
        downsample_docker = downsample_docker,
        start_depth = start_depth,
        final_depth = final_depth_2x,
        seed_override = seed_override,
        reference_fasta = reference_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        original_cram_or_bam_file_read_groups = downSampling_01.read_groups,
        intervals_exons = intervals_exons,
        sample_ID = sample_ID,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        runtime_attr_random_sample = runtime_attr_random_sample,
        runtime_attr_realign = runtime_attr_realign,
        runtime_attr_add_read_group = runtime_attr_add_read_group,
        runtime_attr_mark_duplicates = runtime_attr_mark_duplicates,
        runtime_attr_sort_index = runtime_attr_sort_index,
        runtime_attr_count_coverage = runtime_attr_count_coverage,
        runtime_attr_collect_counts = runtime_attr_collect_counts
    }
  }

  if (run_downsample_4x) {
    call ds2.downSampling_02 as downSampling_02_4x {
      input :
        fastq_file_1 = downSampling_01.fastq_1,
        fastq_file_2 = downSampling_01.fastq_2,
        downsample_docker = downsample_docker,
        start_depth = start_depth,
        final_depth = final_depth_4x,
        seed_override = seed_override,
        reference_fasta = reference_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        original_cram_or_bam_file_read_groups = downSampling_01.read_groups,
        intervals_exons = intervals_exons,
        sample_ID = sample_ID,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        runtime_attr_random_sample = runtime_attr_random_sample,
        runtime_attr_realign = runtime_attr_realign,
        runtime_attr_add_read_group = runtime_attr_add_read_group,
        runtime_attr_mark_duplicates = runtime_attr_mark_duplicates,
        runtime_attr_sort_index = runtime_attr_sort_index,
        runtime_attr_count_coverage = runtime_attr_count_coverage,
        runtime_attr_collect_counts = runtime_attr_collect_counts
    }
  }

  if (run_downsample_6x) {
    call ds2.downSampling_02 as downSampling_02_6x {
      input :
        fastq_file_1 = downSampling_01.fastq_1,
        fastq_file_2 = downSampling_01.fastq_2,
        downsample_docker = downsample_docker,
        start_depth = start_depth,
        final_depth = final_depth_6x,
        seed_override = seed_override,
        reference_fasta = reference_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        original_cram_or_bam_file_read_groups = downSampling_01.read_groups,
        intervals_exons = intervals_exons,
        sample_ID = sample_ID,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        runtime_attr_random_sample = runtime_attr_random_sample,
        runtime_attr_realign = runtime_attr_realign,
        runtime_attr_add_read_group = runtime_attr_add_read_group,
        runtime_attr_mark_duplicates = runtime_attr_mark_duplicates,
        runtime_attr_sort_index = runtime_attr_sort_index,
        runtime_attr_count_coverage = runtime_attr_count_coverage,
        runtime_attr_collect_counts = runtime_attr_collect_counts

    }
  }

  if (run_downsample_8x) {
    call ds2.downSampling_02 as downSampling_02_8x {
      input :
        fastq_file_1 = downSampling_01.fastq_1,
        fastq_file_2 = downSampling_01.fastq_2,
        downsample_docker = downsample_docker,
        start_depth = start_depth,
        final_depth = final_depth_8x,
        seed_override = seed_override,
        reference_fasta = reference_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        original_cram_or_bam_file_read_groups = downSampling_01.read_groups,
        intervals_exons = intervals_exons,
        sample_ID = sample_ID,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        runtime_attr_random_sample = runtime_attr_random_sample,
        runtime_attr_realign = runtime_attr_realign,
        runtime_attr_add_read_group = runtime_attr_add_read_group,
        runtime_attr_mark_duplicates = runtime_attr_mark_duplicates,
        runtime_attr_sort_index = runtime_attr_sort_index,
        runtime_attr_count_coverage = runtime_attr_count_coverage,
        runtime_attr_collect_counts = runtime_attr_collect_counts
    }
  }

  if (run_downsample_custom && defined(final_depth_custom)) {
    call ds2.downSampling_02 as downSampling_02_custom {
      input :
        fastq_file_1 = downSampling_01.fastq_1,
        fastq_file_2 = downSampling_01.fastq_2,
        downsample_docker = downsample_docker,
        start_depth = start_depth,
        final_depth = final_depth_custom,
        seed_override = seed_override,
        reference_fasta = reference_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        original_cram_or_bam_file_read_groups = downSampling_01.read_groups,
        intervals_exons = intervals_exons,
        sample_ID = sample_ID,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        runtime_attr_random_sample = runtime_attr_random_sample,
        runtime_attr_realign = runtime_attr_realign,
        runtime_attr_add_read_group = runtime_attr_add_read_group,
        runtime_attr_mark_duplicates = runtime_attr_mark_duplicates,
        runtime_attr_sort_index = runtime_attr_sort_index,
        runtime_attr_count_coverage = runtime_attr_count_coverage,
        runtime_attr_collect_counts = runtime_attr_collect_counts
    }
  }

  output {
    File fastq_1 = downSampling_01.fastq_1
    File fastq_2 = downSampling_01.fastq_2
    File read_groups = downSampling_01.read_groups

    File? markdup_metrics_2x = downSampling_02_2x.markdup_metrics
    File? sorted_cram_2x = downSampling_02_2x.sorted_cram
    File? crai_file_2x = downSampling_02_2x.crai_file
    File? wgs_coverage_metrics_2x = downSampling_02_2x.wgs_coverage_metrics
    File? read_counts_2x = downSampling_02_2x.read_counts

    File? markdup_metrics_4x = downSampling_02_4x.markdup_metrics
    File? sorted_cram_4x = downSampling_02_4x.sorted_cram
    File? crai_file_4x = downSampling_02_4x.crai_file
    File? wgs_coverage_metrics_4x = downSampling_02_4x.wgs_coverage_metrics
    File? read_counts_4x = downSampling_02_4x.read_counts

    File? markdup_metrics_6x = downSampling_02_6x.markdup_metrics
    File? sorted_cram_6x = downSampling_02_6x.sorted_cram
    File? crai_file_6x = downSampling_02_6x.crai_file
    File? wgs_coverage_metrics_6x = downSampling_02_6x.wgs_coverage_metrics
    File? read_counts_6x = downSampling_02_6x.read_counts

    File? markdup_metrics_8x = downSampling_02_8x.markdup_metrics
    File? sorted_cram_8x = downSampling_02_8x.sorted_cram
    File? crai_file_8x = downSampling_02_8x.crai_file
    File? wgs_coverage_metrics_8x = downSampling_02_8x.wgs_coverage_metrics
    File? read_counts_8x = downSampling_02_8x.read_counts

    File? markdup_metrics_custom = downSampling_02_custom.markdup_metrics
    File? sorted_cram_custom = downSampling_02_custom.sorted_cram
    File? crai_file_custom = downSampling_02_custom.crai_file
    File? wgs_coverage_metrics_custom = downSampling_02_custom.wgs_coverage_metrics
    File? read_counts_custom = downSampling_02_custom.read_counts
  }

}



