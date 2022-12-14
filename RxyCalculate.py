"""https://www.diva-portal.org/smash/get/diva2:1622230/FULLTEXT01.pdf"""
import dadi, sys, os, re, pickle, gzip
from multiprocessing import Pool
import numpy as np
import pylab, sys
import matplotlib.pyplot as plt
import dadi.DFE as DFE
from collections import defaultdict


### let's define a few things: rootdir is where our results files will be created
rootdir = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/'
### Calculate Rxy for the provided input species. Input species is also how directories are organised 
species = sys.argv[1]
sp = species[0]+species.split('_')[1]
speciesdir = '/Users/jennyjames/Dropbox/LASCOUX/Species/'

dir = rootdir+species+'/'

### Specify input vcfs
vcffile_4fold = [file for file in os.listdir(dir) if '_4fold_p50.polarized.vcf.gz' in file][0]
vcffile_0fold = [file for file in os.listdir(dir) if '_0fold_p50.polarized.vcf.gz' in file][0]

print(vcffile_4fold)

### Define populations per species, using a dadi-style popfile found at:
popfile = speciesdir+sp+'/'+sp+'_pops.txt'
print(popfile)
### which population will be used as the reference, to compare the focal populations to?
Refpop = 'FR'

with open(popfile, 'r') as pop_file:
	pop_file = pop_file.readlines()
	pops = set([line.split('\t')[1].split('\n')[0] for line in pop_file])
	pops = [pop for pop in pops if pop != 'GB']
	print(pops)
	

##make VCF data dicts, one entry per site
dd = dadi.Misc.make_data_dict_vcf(dir+vcffile_4fold, popfile)
dd_0 = dadi.Misc.make_data_dict_vcf(dir+vcffile_0fold, popfile)



"""If we were just to look at unique mutations, we could use this code"""
# four_fs = dadi.Spectrum.from_data_dict(dd, ['GB', pop], projections = [40, 40], polarized = True)
# zero_fs = dadi.Spectrum.from_data_dict(dd_0, ['GB', pop], projections = [40, 40], polarized = True)
# 
# # print(four_fs)
# fourfoldFRcounts = four_fs[0:40, 0]
# fourfoldpopcounts = four_fs[0,0:40]
# 
# zerofoldFRcounts = zero_fs[0:40, 0]
# zerofoldpopcounts = zero_fs[0,0:40]
# 
# ### account for frequency of mutant allele- just considering mutations only present in one or the other population
# def Rx(popcounts):
# 	count = 0
# 	for x in list(popcounts):
# 		freq = list(popcounts).index(x)
# 		if x != '--':
# 			ratio_val = list(popcounts).index(x)/len(popcounts)
# 			count = count + (ratio_val * x)
# 	return count
# 
# print(sp + '\t' + pop + '\t' + str(Rx(fourfoldFRcounts)) + '\t' + str(Rx(fourfoldpopcounts)) + '\t' + str(Rx(fourfoldpopcounts)/Rx(fourfoldFRcounts)) + '\t' + str(Rx(zerofoldFRcounts)) + '\t' + str(Rx(zerofoldpopcounts)) + '\t' + str(Rx(zerofoldpopcounts)/Rx(zerofoldFRcounts)) ) 
#print(sp + '\t' + pop + '\t' + str(np.sum(fourfoldFRcounts)) + '\t' + str(np.sum(fourfoldpopcounts)) + '\t'+ str(np.sum(fourfoldFRcounts)/np.sum(fourfoldpopcounts)) + '\t' + str(np.sum(zerofoldFRcounts)) + '\t' + str(np.sum(zerofoldpopcounts)) + '\t' + str(np.sum(zerofoldFRcounts)/np.sum(zerofoldpopcounts)))


### < 1- mutations are efficiently purged relative to FR, > 1 mutations are inefficiently purged relative to FR.

"""Formula requires us to take into account mutant sites that are not strictly unique to either population"""
def freq_pop_calculator(vcfdict, pop):
	freqpopX = 0
	freqpopY = 0
	for k, v in vcfdict.items():
		###are the mutations derived?
		if v['calls'][Refpop][1] != 0 and v['calls'][pop][1] != 0:
	
			###return the mutation counts for our population pairs of interest: format is calls_dict[pop] = (refcalls, altcalls)
			refX, altX = v['calls'][pop]
			refY, altY = v['calls'][Refpop]
			nX = float(refX) + float(altX)
			nY = float(refY) + float(altY)
			fX = float(refX)/nX
			fY = float(refY)/nY
			freqpopX = freqpopX + fX*(1-fY)
			freqpopY = freqpopY + fY*(1-fX)

	return freqpopX, freqpopY

	
for pop in pops:
	PopX, PopY = freq_pop_calculator(dd, pop)
	PopX0, PopY0 = freq_pop_calculator(dd_0, pop)

	print(PopX/PopY)
	print(PopX0/PopY0)


	### new bootstrap- chunk by genomic region, 2Mb
	###however, if there are lots of small scaffolds, this generates lots of small chunks
	### we need to combine them to make 100 chunks for jackknifing- first using partition,
	### then chunker
	Nboot, chunk_size = 100, 2e6
	chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
	chunks0 = dadi.Misc.fragment_data_dict(dd_0, chunk_size)
	
	def partition(lst, n): 
		division = len(lst) / float(n) 
		return [ lst[int(round(division * i)): int(round(division * (i + 1)))] for i in range(n) ]

	###list of 100 dictionaries- unpack into single dictionary
	chunks = partition(chunks, Nboot)
	chunks0 = partition(chunks0, Nboot)

	def dict_chunker(list_of_dicts):
		dd_chunks = []
		for chunk in list_of_dicts:			
			dd_100_chunks = {}
			for key, value in chunk[0].items():		
				dd_100_chunks[key] = value		
			dd_chunks.append(dd_100_chunks)
		return(dd_chunks)
	
	dd_chunks = dict_chunker(chunks)
	dd_0_chunks = dict_chunker(chunks0)
		
		
	with open(rootdir+'Rxy_to_'+Refpop+'_norm/Rxy'+sp+pop+'.tsv', 'a') as resfile:
		resfile.write(str(PopX/PopY) + '\t'+ str(PopY/PopX) + '\t'+ str(PopX0/PopY0) + '\t'+ str(PopY0/PopX0) + '\t' + str((PopX0/PopX)/(PopY0/PopY)) + '\t' + str((PopY0/PopY)/(PopX0/PopX)) +'\n')

		### we will now jackknife, recalculating the Rxy statistic 100 times, excluding 1 chunk each time
		for i in range(0, Nboot):
# 			print(len(dd_chunks))	
# 			print(i)	
			dd_chunks_new = dd_chunks[:]
			dd_0_chunks_new = dd_0_chunks[:]			
			del dd_chunks_new[i]
			del dd_0_chunks_new[i]
			 			
			dd_del = {}
			for chunk in dd_chunks_new:
				for key, value in chunk.items():		
					dd_del[key] = value	

			dd_0_del = {}
			for chunk in dd_0_chunks_new:
				for key, value in chunk.items():		
					dd_0_del[key] = value						
					
			chunkX, chunkY = freq_pop_calculator(dd_del, pop)	
			chunkX0, chunkY0 = freq_pop_calculator(dd_0_del, pop)	
			resfile.write(str(i)+ '\t' + str(chunkX) + '\t' + str(chunkY) + '\t' + str(chunkX/chunkY) + '\t' + str(chunkX0) + '\t' + str(chunkY0) + '\t' + str(chunkX0/chunkY0) + '\t' + str((chunkX0/chunkX)/(chunkY0/chunkY))+ '\t' + str((chunkY0/chunkY)/(chunkX0/chunkX)) + '\n')




























