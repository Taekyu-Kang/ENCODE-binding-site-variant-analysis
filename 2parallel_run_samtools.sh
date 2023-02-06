#!/bin/bash

#runs run_samtools.py for n parallel jobs
#n=3
DIR=$1 #parent directory containing BAM files
AFILES=$2 #metadata on files used
CFILES=$3 #list of completed files, can be NA

for ((i=1; i<=3; i++))

do
	python2 /home/tkang/design/code/run_samtools.py $DIR $AFILES $CFILES $i &
done

