#!/bin/bash
set -euo pipefail
if [ $# -lt 1 ] ; then
    echo ""
    echo "usage: count_fastq.sh [fastq_file1] <fastq_file2> ..|| *.fastq"
    echo "counts the number of reads in a fastq file"
    echo ""
    exit 0
fi

filear=${@};
for i in ${filear[@]}
do
lines=$(wc -l $i|cut -d " " -f 1)
count=$(($lines / 4))
echo "$count"
done


#simple script to count number of reads in a fastq file