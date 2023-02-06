#identifies binding sites located -3kb and +1kb of transcriptional start sites
#outputs a text file containing the following information
#Gene=gene symbol
#TranscriptID=Ensembl transcript ID
#ChIP=ChIP-seq protein used
#CellType=cell type
#BindingSite=ChIP-seq protein binding site
#BSseq=binding site sequence
#SNP=location of SNP, reference/alternative
#BioSample=additional information about the experiment

from sys import argv
import pickle
from Bio import SeqIO

#extracts information from a GTF
def parse_annot(filename,sep="\t"):
	annot = []
	with open(filename) as f:
		for line in f:
			linedata = line.strip().split(sep)
			annot.append(linedata)
	return annot 

#load TF binding site data pickle dict
def load_pkl(filename):
	with open(filename,"rb") as f:
		return pickle.load(f)

#extracts desired information from the metadata file
def parse_file(filename,sep="\t"):
	fdat = {}
	with open(filename) as f:
#		f.readline()
		for line in f:
			linedata = line.strip().split(sep)
			exp = linedata[0]
			tf = linedata[1]
			cell = linedata[2]
			bio = linedata[3]

			fdat[exp+"."+tf] = [cell,bio]
	return fdat

#identifies binding sites located near TSSs
def get_BS(chrdict,annot,fdat):
	final = []
	for i in annot:
		chrom = i[0]
		#generate start and end sites for the location around the TSS
		if i[4] == "+":
			tss = int(i[2])
			start = tss-3000
			end = tss+1000
		elif i[4] == "-":
			tss = int(i[3])
			end = tss+3000
			start = tss-1000
		tID,gene = i[5:]

		if gene == "havana":
			continue

		try: #skips chromosomes not in my genome file
			chrdat = chrdict[chrom]
		except:
			continue

		sw = False
		ew = False
		#I want to make sure I capture as many potential binding sites as possible
		#as I am comparing the window around the TSS to the window keys, I want to make sure that if the TSS window overlaps two window keys, I'll search both for binding sites
		for w in chrdat:
			if w[1] < start:
				continue
			elif start > w[0] and start < w[1]:
				sw = chrdat[w]
				if end > w[0] and end < w[1]:
					ew = sw
					break
				else:
					continue
			elif end > w[0] and end < w[1]:
				ew = chrdat[w]
				break
			else:
				continue

		#sometimes there are no matching windows for the start and/or end
		if not sw and not ew:
			continue
		if not sw:
			sw = []
		if not ew:
			ew = []

		if sw == ew:
			for site in sw:
				if site[0] > start and site[1] < end:
					final.append([tID,gene,tss,chrom,site])
		elif sw != ew:
			for site in sw:
				if site[0] > start and site[1] < end:
					final.append([tID,gene,tss,chrom,site])			
			for site in ew:
				if site[0] > start and site[1] < end:
					final.append([tID,gene,tss,chrom,site])

	return final

#reads in concatenated genome file and returns a dictionary with chromosome as keys and sequences as values
def make_genome_dict(filename):
	genome = {}
	fa_dat = SeqIO.parse(filename,"fasta")
	print "Starting to make genome dict..."

	for chrom in fa_dat:
		genome[chrom.id] = chrom
	print "Done making genome dict,", len(genome), "chromosomes."
	return genome

#takes the output from get_BS and annotates with metadata
#writes output
def annotate(final,fdat,genome):

	wf = open(outdir+"Design_ENCODE_BS_var.txt","w")
	wf.writelines("#Gene=gene symbol\n")
	wf.writelines("#TranscriptID=Ensembl transcript ID\n")
	wf.writelines("#ChIP=ChIP-seq protein used\n")
	wf.writelines("#CellType=cell type\n")
	wf.writelines("#BindingSite=ChIP-seq protein binding site\n")
	wf.writelines("#BSseq=binding site sequence\n")
	wf.writelines("#SNP=location of SNP, reference/alternative\n")
	wf.writelines("#BioSample=additional information about the experiment\n")
	wf.writelines("Gene\tTranscriptID\tChIP\tCellType\tBindingSite\tBSseq\tDist\tSNP\tBioSample\n")

	for i in final:
		start,end,tfID,pos,var = i[4]
		snp = str(pos)+":"+var
		exp,tf = tfID.split(".")
		cell,biosample = fdat[tfID]
		chrom = i[3]
		tss = i[2]
		trID = i[0]
		gene = i[1]
		seq = str(genome[chrom][start:end+1].seq)
		dist = int(start)-int(tss)
		bs = chrom+":"+str(start)+"-"+str(end)

		wf.writelines(gene+"\t"+trID+"\t"+tf+"\t"+cell+"\t"+bs+"\t"+seq+"\t"+str(dist)+"\t"+snp+"\t"+biosample+"\n")

	wf.close()

if __name__=="__main__":
	outdir = argv[1] #desired output directory
	annot_file = argv[2] #GTF
	site_file = argv[3] #BS pickle dict
	genome_file = argv[4] #concatenated genome file
	meta_file = argv[5] #metadata file

	print "Loading files..."
	annot_dat = parse_annot(annot_file)
	site_pkl = load_pkl(site_file)
#	for i in site_pkl:
#		print i, len(site_pkl[i])
	meta_dat = parse_file(meta_file)
	print "Done loading files, starting on binding site analysis."

	BSfinal = get_BS(site_pkl,annot_dat,meta_dat)

	genome_dat = make_genome_dict(genome_file)
	annotate(BSfinal,meta_dat,genome_dat)

	print "Done!"

