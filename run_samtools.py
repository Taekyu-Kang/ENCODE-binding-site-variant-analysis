#runs samtools to get only reads that map uniquely and sort
#can also run depth if possible, couldn't in my case due to space
#called by 2parallel_run_samtools.sh

from sys import argv
import subprocess as sp
import os

#first get metadata of files used, split into n lists for parallel run
#returns a dictionary with a number 1 through n as a key and lists as values
def parse_file(filename,sep="\t"):
	files = []
	splitn = int(keyn)
	with open(filename) as f:
#		f.readline()
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

#	for i in fdict:
#		print i, len(fdict[i])

	return fdict

#follow up for if the script stalls or crashes, so I can skip the files that are already completed
def parse_complete(filename,sep="\t"):
	complete = []
	with open(filename) as f:
		for line in f:
			linedata = line.strip()
			complete.append(linedata)
	return complete

#run samtools to end up with sorted BAM
def run_samtools(fdict,keyn,complete):
	indir1 = indir+"BAM/"
	#grab one of the lists as told by bash file to loop through
	flist = fdict[keyn]
	for i in flist:
		if i in complete:
			continue
		print "Working on", i, "..."
		cmd = ["mkdir", indir+"final/"+i]
		sp.call(cmd)

		outdir = indir+"final/"+i+"/"
		#print outdir

		mapped = sp.Popen(["samtools", "view", "-u", "-F","4", indir1+i+".bam"],stdout=sp.PIPE)
		
		unique = sp.Popen(["samtools", "view", "-u", "-q", "30"],stdin=mapped.stdout,stdout=sp.PIPE)
		
		mapped.stdout.close()
		
		sp.call(["samtools", "sort", "-O", "BAM", "-T", outdir+i+"_temp", "-o", outdir+i+"_sorted.bam"],stdint=unique.stdout)
		
		unique.stdout.close()
		
#		cmd = ["samtools", "depth", "-o",  outdir+i+".depth", outdir+i+"_sorted.bam"]
#		sp.call(cmd)

		final = os.listdir(outdir)
		#remove original BAM file to save space
		if i+"_sorted.bam" in final:
			cmd = ["rm", indir1+i+".bam"]
			sp.call(cmd)

		print "Done with", i, "!"

if __name__=="__main__":
	indir = argv[1] #parent directory containing BAM files
	all_files = argv[2] #metadata on files used
	completed_files = argv[3] #list of completed files, can be NA
	keyn = str(argv[4]) #key to use from bash file

	if completed_files == "NA":
		complete_dat = []

	f_dict = parse_file(all_files)
	run_samtools(f_dict,keyn,complete_dat)
