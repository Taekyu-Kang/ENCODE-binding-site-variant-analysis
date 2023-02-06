#!/bin/bash

#runs run_bcftools.py for n parallel jobs
#n=3
DIR=$1 #parent directory containing BAM files
AFILES=$2 #metadata on files used

for ((i=1; i<=3; i++))

do
	python2 /home/tkang/design/code/run_bcftools.py $DIR $AFILES $i &
done

