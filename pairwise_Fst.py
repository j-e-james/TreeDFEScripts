###  write up the commands to conduct polyDFE independent inference of the DFE on every pairwise comparison
import dadi, sys, os, re, pickle
import numpy as np
import pylab
import matplotlib.pyplot as plt
import dadi.DFE as DFE
from collections import defaultdict

popfile = '/Users/jennyjames/Dropbox/LASCOUX/Species/Pabies/Pabies_pops.txt'
# annofile = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/gentree_annotation/PS.v3.full.annotation.allSites.txt'
vcfdir = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/Picea_abies/'

sample_size = 40

vcffile_4fold = [x for x in os.listdir(vcfdir) if '4fold_p50.polarized.vcf.gz' in x]
vcffile_0fold = [x for x in os.listdir(vcfdir) if '0fold_p50.polarized.vcf.gz' in x]
vcffile_4fold = vcfdir+vcffile_4fold[0]
vcffile_0fold = vcfdir+vcffile_0fold[0]
			

workdir = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/Picea_abies'

filelist = []

for file in os.listdir(workdir):
	if '_40.in' in file and '.out' not in file:
		print(file)
		filelist.append(workdir+file)

dd = dadi.Misc.make_data_dict_vcf(vcffile_4fold, popfile)
dd_0 = dadi.Misc.make_data_dict_vcf(vcffile_0fold, popfile)

for focus_file in filelist:
	popfilelist = [x for x in filelist if x != focus_file]
	for file in popfilelist:
		focus_file_pop = focus_file[-8:][:2]
		file_pop = file[-8:][:2]
# 		print(focus_file_pop+' '+file_pop)
		fs = dadi.Spectrum.from_data_dict(dd, [focus_file_pop, file_pop], projections = [sample_size, sample_size], polarized = True)
		fs_0 = dadi.Spectrum.from_data_dict(dd_0, [focus_file_pop, file_pop], projections = [sample_size, sample_size], polarized = True)
		print("c('Picea_abies','"+focus_file_pop+"-"+file_pop+"','"+focus_file_pop+"',"+str(fs.Fst())+","+str(fs_0.Fst())+"),")

	

