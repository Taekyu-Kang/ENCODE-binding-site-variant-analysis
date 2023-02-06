ENCODE binding site variant analysis README

This pipeline is used to process data from the ENCODE database to identify binding sites located around gene TSSs containing SNPs against the reference. 

Script 1 downloads all the necessary files and organizes them by experiment. Returns a metadata sheet of all experiments used.
Scripts 2-3 processes BAM files using samtools and bcftools in parallel, calls run_samtools.py and run_bcftools.py respectively. Returns VCFs of filtered SNPs per experiment.
Script 4 combines VCF and BED files and filters for binding sites with a called SNP.
Script 5 combines all files for all experiments, returns a searchable dictionary with genomic locations as keys and binding sites as values.
Script 6 takes output from Script 5, combines with annotations to deliver a text file containing binding sites containing SNPs located around TSSs.
