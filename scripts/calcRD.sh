#!/bin/bash
err=$'No arguments supplied. The wanted down sampling depths required:\nArg1 = Original Read Depth\nArg2 = Desired Final Read Depth\nArg3 = Original Read Counts'

if [ $# -eq 0 ]
  then
    echo "$err"
fi

ord=$1
frd=$2
count=$3

numerator=$(bc <<< "$frd*$count")
#numerator=$(($count * $frd))
#echo $numerator

readsToSample=$(bc <<< "$numerator / $ord")

echo $readsToSample