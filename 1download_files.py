#script to download bam and bed files for the experiments in the metadata file
#also generates new metadata file of files used

from sys import argv
import os
import subprocess as sp

#extract information from the metadata file from ENCODE cart
def parse_meta(filename,sep="\t"):
	all_files = {}
	with open(filename) as f:
		f.readline()
		for line in f:
			if line.startswith("ID"):
				linedata = line.strip().split(sep)
				kl = linedata
			else:
				linedata = line.strip().split(sep)
				if len(linedata) < len(kl):
					for i in range(len(kl)-len(linedata)):
						linedata.append("")
				#print linedata, len(linedata)
				temp_dict = {kl[i]:linedata[i] for i in range(len(kl))}
				#print temp_dict
				
				#get simplified version of Biosample column
				all_files[temp_dict["Accession"]] = [temp_dict["Target of assay"],temp_dict["Biosample term name"], ";".join([i+":"+temp_dict[i] for i in temp_dict if i in ["Organism","Life stage","Biosample age","Biosample treatment"]])]
				t_files = [i.split("/")[2] for i in temp_dict["Files"].split(",")]
				all_files[temp_dict["Accession"]].append(t_files)

	print len(all_files), "files from ENCODE search."
	return all_files

#in case the script crashes, get list of completed files
def parse_done(filename,sep="\t"):
	done = []
	with open(filename) as f:
		f.readline()
		for line in f:
			linedata = line.strip().split(sep)
			done.append(linedata[0])
	print len(done), "files already done."
	return done

#in the metadata file, I wasn't sure which file ID was the BED or BAM files, so I looped through them until I got the right file
def download(all_files,done):

	if stopped == "No":
		wf = open(outdir+"ENCODE_TF_files_used.txt","w")
		wf.writelines("Accession\tTarget_of_assay\tBiosample_term_name\tInfo\tBED_file\tBAM_file\n")
		wf.close()

	print "Starting file extraction from ENCODE..."
	for x in all_files:
		if x in done:
			continue
		xdat0 = all_files[x]
		bed = False
		bam = False
		print x
		xdat = xdat0[3]
		for f in xdat:
			if not bed:
				bed_cmd = ["wget","https://www.encodeproject.org/files/"+f+"/@@download/"+f+".bed.gz"]
				bed_out = sp.call(bed_cmd)
				#print bed_out
				if bed_out == 0:
					bed = f
			if not bam:
				bam_cmd = ["wget","https://www.encodeproject.org/files/"+f+"/@@download/"+f+".bam"]
				bam_out = sp.call(bam_cmd)
				#print bam_out
				if bam_out == 0:
					bam = f

		print "****"
		print "****"
		print x, bed, bam

		if not bed or not bam:
			continue

		with open(outdir+"ENCODE_TF_files_used.txt","a") as wf:
			wf.writelines(x+"\t"+"\t".join(all_files[x][:3])+"\t"+bed+"\t"+bam+"\n")


if __name__=="__main__":
	outdir = argv[1] #desired output directory
	metafile = argv[2] #cart output from ENCODE
	stopped = argv[3] #check if stopped, Yes or No

	if stopped == "Yes":
		done_dat = parse_done(outdir+"ENCODE_TF_files_used.txt")
	elif stopped == "No":
		done_dat = []

	meta_dat = parse_meta(metafile)
	download(meta_dat,done_dat)

	print "Done!"


