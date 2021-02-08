#!/bin/bash
set -euo pipefail
err=$'No arguments supplied. The wanted down sampling depths required:\nArg1 = Original Read Depth\nArg2 = Desired Final Read Depth\nArg3 = Original Read Counts'

if [ $# -eq 0 ]
  then
    echo "$err"
fi

ord=$1
frd=$2
count=$3

numerator=$(echo $frd $count| awk ' { printf "%0.2f\n", ($1 * $2); } ')

readsToSample=$(echo $numerator $ord| awk ' { printf "%0.2f\n", ($1 / $2); } ')

echo $readsToSample