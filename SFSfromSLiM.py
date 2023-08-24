import re, numpy, os, pandas


"""
Next is the Mutations section, which lists all of the currently segregating mutations in the
population. Each mutation is listed on a separate line. The first field is a within-file numeric
identifier for the mutation, beginning at 0 and counting up (although mutations are not listed in
sorted order according to this value); see below for a note on why this field exists. Second is the
mutation’s id property (see section 24.10.1), a within-run unique identifier for mutations that does
not change over time, and can thus be used to match up information on the same mutation within
multiple output dumps made at different times. Third is the identifier of the mutation’s mutation
type, such as m1. Fourth is the position of the mutation on the chromosome, as a zero-based base
position. Fifth is the selection coefficient of the mutation, and sixth is its dominance coefficient
(the latter being a property of the mutation type, in fact). Seventh is the identifier for the
subpopulation in which the mutation originated, and eighth is the tick in which it originated.
Finally, the ninth field gives the mutation’s prevalence: an integer count of the number of times that
it occurs in any genome in the population. 
"""


def SFS_convert(pop_mutations, n):
	###converts SLiM mutation output into SFS
	SFS = {} # full SFS with empty sites
	for x in range(1, n+1):
		SFS[x] = []
	for mut in pop_mutations:
		SFS[int(mut)].append(mut)
	SFS_counts = []
	for k, v in SFS.items():
		SFS_counts.append(len(v))
	return SFS_counts

SFS_neutral = []
SFS_deleterious = []
s_values = [] 
for file in os.listdir('/Users/jennyjames/Desktop/GeneralScripts'):
	if '-2_' in file:

		with open(file, 'r') as SLiMOutput:
			SLiMOutput = SLiMOutput.read()

			sample_number = int(SLiMOutput.split('#OUT: ')[1].split('\n')[0].split(' ')[4])
# 			print(sample_number)

			Ratio_syn_nonsyn = SLiMOutput.split('initializeGenomicElementType(')[1].split(';')[0]
			Ratio_syn_nonsyn = [x.split(',') for x in re.findall(r'\((.*?)\)', Ratio_syn_nonsyn)]
			Mutation_types, Ratios = Ratio_syn_nonsyn[0], Ratio_syn_nonsyn[1] 	
			Ratios = [int(x) for x in Ratios]
			Ratios = [x/numpy.sum(Ratios) for x in Ratios]

			Length = SLiMOutput.split('initializeGenomicElement(')[1].split(');')[0]
			Length = int(Length.split(', ')[-1])+1
			Lengths = [Length*x for x in Ratios]
	
# 			print(Mutation_types)
# 			print(Lengths)
	
			Mut_freq_list = []
	
			Mutations = SLiMOutput.split('Mutations:')[1].split('Genomes:')[0].split('\n')
			Mutations = [x for x in Mutations if len(x) != 0]
			
			for x in Mutation_types:
				x_Mutations = [y for y in Mutations if x in y]
				if x != 'm1':
					s_values.append(', '.join([str(z.split(' ')[4]) for z in x_Mutations]))
					
				
				Mut_freq_list.append([z.split(' ')[-1] for z in x_Mutations])
			
			
			SFS_neutral.append(SFS_convert(Mut_freq_list[0], sample_number))
			SFS_deleterious.append(SFS_convert(Mut_freq_list[1], sample_number))
			
# print('\n'.join(SFS_neutral))
# print('\n\n\n')
# print(SFS_deleterious)

SFS_neutral_sum = []
for x in zip(*SFS_neutral):
	SFS_neutral_sum.append(str(numpy.mean(x)))
	
SFS_deleterious_sum = []
for x in zip(*SFS_deleterious):
	SFS_deleterious_sum.append(str(numpy.mean(x)))

print('SFS <- data.frame(c('+', '.join(str(x) for x in range(1, 500+1))+'),')
print('c('+', '.join(SFS_neutral_sum)+'),')	
print('c('+', '.join(SFS_deleterious_sum)+'))')	
print('s_values <- data.frame(c('+', '.join(s_values)+'\n))')




# print('\n\n\n')
# print(', '.join(s_values))

						
			### for SFS:
# 			for mutation_set, Ls in zip(Mut_freq_list, Lengths):
# 				print(', '.join(SFS_convert(mutation_set, sample_number))
	
	
	

	
	