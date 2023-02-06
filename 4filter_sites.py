#combines BED and VCF for each experiment and fitler for binding sites with a variant

from sys import argv
import os
import subprocess as sp

#metadata on files used
#this time returns dictionary with experiment ID as the key with files used and TF name as values
def parse_file(filename,sep="\t"):
	files = {}
	with open(filename) as f:
		f.readline()
		for line in f:
			linedata = line.strip().split(sep)
			bam = linedata[-1]
			bed = linedata[-2]
			exp = linedata[0]
			tf = linedata[1]

			files[exp] = [bam,bed,tf]

	return files

#run bedtools intersect to get binding sites containing variants
def filter_sites(files):
	for i in files:
		fdat = files[i]
		bam, bed, tf = fdat

		cmd = ["bedtools", "intersect", "-a", indir+"BED/"+bed+".bed", 
		"-b", indir+"final/"+bam+"/"+bam+".filtsnps.vcf.gz", "-wa", "-wb"]
		wf = open(indir+"filtered/"+i+"."+tf+".txt","w")
		sp.call(cmd,stdout=wf)
		wf.close()


if __name__=="__main__":
	indir = argv[1] #parent directory containing output
	all_files = argv[2] #metadata on files used

	f_dict = parse_file(all_files)
	filter_sites(f_dict)
