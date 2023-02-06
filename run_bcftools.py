#run bcftools mpileup and call to generate VCF, then view to filter for SNPs
#called by 3parallel_run_bcftools.sh

from sys import argv
import subprocess as sp
import os

#first get metadata of files used, split into n lists for parallel run
#returns a dictionary with a number 1 through n as a key and lists as values
def parse_file(filename,sep="\t"):
	files = []
	splitn = 3
	with open(filename) as f:
		f.readline()
		for line in f:
			linedata = line.strip().split(sep)
			bam = linedata[-1]
			files.append(bam)

	fdict = {}
	ct = 0
	splits = len(files)/splitn
	for i in range(0,len(files),splits):
		ct += 1
		fdict[str(ct)] = files[i:i+splits]

	return fdict

#run bcftools to get SNPs
def run_bcftools(fdict,keyn):
	#grab one of the lists as told by bash file to loop through
	flist = fdict[keyn]
	for i in flist:

		print "Working on", i, "..."

		outdir = indir+"final/"+i+"/"
		#print outdir

		mpileup = sp.Popen(["bcftools", "mpileup", "-Ou", "-f", "/home/tkang/design/hg/Homo_sapiens.GRCh38.dna.cat.fa", outdir+i+"_sorted.bam"], stdout=sp.PIPE)

		bcall = sp.Popen(["bcftools", "call", "-mv", "-Ou"], stdin=mpileup.stdout,stdout=sp.PIPE)
		
		mpileup.stdout.close()
		
		sp.call(["bcftools", "view", "-I", "-V", "indels", 
		"-e", "QUAL<50", "-Oz", 
		"-o", outdir+i+".filtsnps.vcf.gz"],stdin=bcall.stdout)
		
		bcall.stdout.close()
		
		final = os.listdir(outdir)
		#remove sorted file to save space
		if i+".vcf.gz" in final:
			cmd = ["rm", outdir+i+"_sorted.bam"]
			sp.call(cmd)

		print "Done with", i, "!"

if __name__=="__main__":
	indir = argv[1] #parent directory containing BAM files
	all_files = argv[2] #metadata on files used
	keyn = str(argv[3]) #key to use from bash file

	f_dict = parse_file(all_files)
	run_bcftools(f_dict,keyn)
