import gzip, os, re, dadi, numpy, sys
from collections import defaultdict
from collections import Counter

	
def polyDFESFS(line):
	fs = line.split('\t')[:-1]
	SFS = {}
	for x, y in enumerate(fs):
		if y != 'masked':
			SFS[x] = y
	n = len(fs)+1
	return SFS, n

	
### caclulate pi, given an SFS of format freq, no. mutations of that freq, a Len (number of sites), and the number of genomes
def PiCalculate(SFS, Len, subpop_genomes):
	count = 0
	for x in SFS:
		count = count + float(SFS[x])*2* (float(x)/float(subpop_genomes)*(1-(float(x)/float(subpop_genomes))))
	pi = count/Len * float(subpop_genomes)/(float(subpop_genomes)-1)
	return pi	


sfs_filepath = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/'
sp = sys.argv[1]
# sfs_file = sfs_filepath+sp+'/'+sys.argv[2]

sfs_file = sys.argv[2]

with open(sfs_file, 'r') as sfs:
	sfs = sfs.readlines()
	SFS4 = sfs[1]
	SFS0 = sfs[2]
	
	zerofoldcount = int(SFS0.split('\t')[-1])
	fourfoldcount = int(SFS4.split('\t')[-1])

	SFS4, genomes4 = polyDFESFS(SFS4)
	SFS0, genomes0 = polyDFESFS(SFS0)
	
	
	pi4fold = PiCalculate(SFS4, fourfoldcount, genomes4)
	pi0fold = PiCalculate(SFS0, zerofoldcount, genomes0)
	print(sys.argv[1]+'\t'+sys.argv[2]+'\t'+str(pi4fold)+'\t'+str(pi0fold)+'\t'+str(pi0fold/pi4fold))
	
	
	
	
	
	
	
	
	
	
	
	
	
	