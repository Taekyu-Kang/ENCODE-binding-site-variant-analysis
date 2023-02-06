#combines all the sites for all experiments used into one file
#outputs a dictionary with windows of genomic locations as keys and binding sites that fall in those windows as values
#reduces the time spent searching through each position when looking for sites located around TSSs

from sys import argv
import collections
import os
import pickle

#had to check the lengths of the rows for each file, sometimes the bed files had different column counts

def check_len(filename,sep="\t"):
	tf = filename.split(".txt")[0]
	with open(indir+filename) as f:
		for line in f:
			linedata = line.strip().split(sep)
			ll = len(linedata)
	print tf, ll, linedata

#extract the binding site and SNP info from each experiment
#outputs to global chr_dict
def parse_file(filename,sep="\t"):
	tf = filename.split(".txt")[0]
	with open(indir+filename) as f:
		for line in f:
			linedata = line.strip().split(sep)
			splitdat = []
			#as row lengths could be different, I had to check when the data for the BED vs VCF started to extract the information properly
			for i in range(len(linedata)):
				if linedata[i].startswith("chr"):
					splitdat.append(linedata[i:i+5])
			#print splitdat
			if linedata[0] not in chr_dict:
				chr_dict[linedata[0]] = []
			chr_dict[linedata[0]].append([int(linedata[1]),int(linedata[2]),
				tf,int(splitdat[1][1]),"/".join(splitdat[1][3:5])])
				
	print "Done with", filename

#sorts all the binding sites for each chromosome, then generates windows to use as keys in the location dictionary
#generating 2000 windows for each chromosome
def make_windows(chrdat):
	tempdict = {}
	all_pos = []
	for i in chrdat:
		all_pos.append(i[0])
		all_pos.append(i[1])
	maxx = max(all_pos)
	minn = min(all_pos)
	splits = (maxx-minn)/2000
	for i in range(minn,maxx,splits):
		tempdict[(i,i+splits)] = []
	tempdict1 = collections.OrderedDict(sorted(tempdict.items()))

	for i in chrdat:
		for w in tempdict1:
			if i[0] < w[0]:
				continue
			elif i[0] > w[0] and i[0] < w[1]:
				tempdict1[w].append(i)
				break

	return tempdict1

if __name__=="__main__":
	indir = argv[1] #parent directory containing output files
	outdir = argv[2] #desired output directory
	all_files = os.listdir(indir)
	chr_dict = {}
	for file in all_files:
		#check_len(file)
		parse_file(file)
	for chrom in chr_dict:
		print chrom, len(chr_dict[chrom])

	print "Starting on making windows..."
	new_dict = {}

	for chrom in chr_dict:
		chrom_dat = chr_dict[chrom]
		new_dat = make_windows(chrom_dat)
		new_dict[chrom] = new_dat
		print "Done with", chrom, "!"

	print new_dict.keys()
	#saves the dictionary as a pickle object
	with open(outdir+"all_TF_withVAR_byCHR.pkl","wb") as f:
		pickle.dump(new_dict,f,pickle.HIGHEST_PROTOCOL)
	print "DONE!"
