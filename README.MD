<<<<<<< HEAD
<h2><b>Genome Downsampling Wokflow</b></h2>
=======
<h2><b>Genome Downsampling Wokflow<b/></h2>
>>>>>>> 3c1dfb85e757b1794876701572084bf5bd798c51

The genome downsampling workflow consists of <b>downsampling_part1</b> and <b>downsampling_part2</b> and executed in Downsample.wdl.

<b>Downsampling_part1:</b>
Converts a CRAM/BAM input file into paired-end FASTQ files and extracts read groups from input file.

<b>Downsampling_part2:</b>
Handles downsampling to a specified final coverage and assumes an initial 30x coverage. It outputs read counts, downsampled CRAM file, coverage estimate, marked duplicate reads.
