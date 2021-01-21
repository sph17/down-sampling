<h2><b>Genome Downsampling Wokflow</b></h2>

The genome downsampling workflow consists of <b>downsampling_part1</b> and <b>downsampling_part2</b> and executed from <b>Downsample.wdl</b>. The workflow utilizes <a href="https://github.com/samtools/samtools">Samtools</a>, <a href="https://github.com/broadinstitute/picard">Picard</a>, <a href="https://sourceforge.net/projects/bio-bwa/files/">BWA</a>, <a href="https://github.com/broadinstitute/gatk/releases/tag/4.1.9.0">GATK</a>, and <a href="https://github.com/dcjones/fastq-tools">fastq-tools</a> to randomly down-sample to a desired sequencing coverage.

<b>Downsampling_part1:</b>
Converts a CRAM/BAM input file into paired-end FASTQ files and extracts read groups from input file.

<b>Downsampling_part2:</b>
Handles downsampling to a specified final coverage and assumes an initial 30x coverage. It outputs read counts, downsampled CRAM file, coverage estimate, marked duplicate reads.

<i>Note:</i> Dockerfile installs: 
<ul>
<li>bwa version 0.7.17, which is different than publically available bwa internal index files if using bwakit (0.7.12).</li>
<li>samtools 1.11, which may require regeneration of faidx</li>
<li>picard 2.23.8</li>
<li>fastq-tools 0.8.3</li>
</ul>
