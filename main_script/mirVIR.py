# easily extendable script that does all of my analyses

import string
import re
import scipy
import time
import math
import random
import operator

import numpy as np
from numpy.random import randn
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


import sys, os, string, numpy, matplotlib.pyplot as plt
from scipy.stats import spearmanr
from distance import hamming
from numpy import mean



def breakfile(line):
	temp = line
	temp = string.replace(temp, "\n", "")
	temp = temp.split("\t")

	return temp


def bincountdict(lst):
	countdict = {}
	for item in lst:
		if item in countdict:
			countdict[item] += 1
		else:
			countdict[item] = 1

	return countdict


def bincount(xlab,name,lst,filename):
	data = lst[:]
	alldata = sorted(data[:])
	labs = sorted(list(set(data[:])))
	howmany = {}
	for el in labs:
		howmany[el] = 0

	for el in alldata:
		howmany[el] += 1





	bins = []

	for el in labs:
		bins.append(howmany[el])

	left = range(len(bins))

	plt.bar(left, bins)
	plt.xticks(left,labs,rotation=35)
	plt.title(name)
	plt.xlabel(xlab)

	plt.savefig(filename, bbox_inches='tight')
	plt.close()
	return


def boxplot(bin_forming_lst,bin_elements, title, x_axis, y_axis, save_name):
	bins_4_pic = {}

	for index, item in enumerate(bin_elements):
		bins_4_pic.setdefault(float(bin_forming_lst[index]), []).append(item)
	labels = sorted(bins_4_pic.keys())
	nums = [bins_4_pic[i] for i in labels]

	plt.figure(figsize=(10,7))
	plt.boxplot(nums,labels=labels)
	plt.xlabel(x_axis)
	plt.ylabel(y_axis)
	plt.title(title)
	plt.savefig("images/"  + save_name + '.png')
	plt.close()


def mirna_rates(mirna2age):
	mirna_ages = sorted(list(set(mirna2age.values())))
	number_in_agecount = {}
	for age in mirna_ages:
		number_in_agecount[age] = 0

	for age in mirna2age.values():
		number_in_agecount[age] += 1

	rates_lst = {}
	fle = open('rates.txt','w')
	for i in range(1,len(mirna_ages)):
		end_age = mirna_ages[i-1]
		start_age = mirna_ages[i]

		calc = (float(number_in_agecount[end_age]) / float(start_age - end_age))
		key = '%f to %f' %(start_age, end_age)
		fle.write(key + ':' + str(calc) + '\n')
		rates_lst[key] = calc

	fle.close()

	fle = open('txtfles/mirnas2age.txt','w')
	
	sorted_ages = sorted(mirna2age.items(), key=operator.itemgetter(1))


	for item in sorted_ages:
		fle.write(str(item[1]) + '\t' + item[0] + '\n')
	fle.close()




	








def break_files(mirna_ages,gene_ages,diseases,family_associations):
	fle1 = open(mirna_ages,'r')
	mirna_ages_unparsed = fle1.readlines()
	fle1.close()

	fle2 = open(gene_ages,'r')
	gene_ages_unparsed = fle2.readlines()
	fle2.close()

	fle3 = open(diseases,'r')
	disease_associations = fle3.readlines()
	fle3.close()

	fle4 = open(family_associations,'r')
	families = fle4.readlines()
	fle4.close()

	gene2age = {}

	for line in gene_ages_unparsed:
		p = breakfile(line)
		gene2age[p[0]] = float(p[1])


	mirna2age = {}
	age2mirna = {}

	for line in mirna_ages_unparsed:
		if line[0] == '#':
			continue
		p = breakfile(line)
		mirna2age[p[0]] = float(p[1])
		age2mirna.setdefault(float(p[1]),[]).append(p[0])

	disease2mirna = {}
	disease2age = {}
	mirna2disease = {}
	age2disease = {}

	for line in disease_associations:
		if line[0] == '#':
			continue
		p = breakfile(line)
		if p[0] not in mirna2age:
			continue
		mirna2disease.setdefault(p[0],[]).append(p[1])
		disease2mirna.setdefault(p[1],[]).append(p[0])

	for disease in disease2mirna:
		mirnalst = disease2mirna[disease]
		disease2age[disease] = [mirna2age[mirna] for mirna in mirnalst]

	for age in age2mirna:
		mirnalst = age2mirna[age]
		diseaselst = []
		for mirna in mirnalst:
			if mirna in mirna2disease:
				diseaselst += mirna2disease[mirna]
		age2disease[age] = list(set(diseaselst))

	family2members = {}
	member2family_name = {}

	for line in families:
		if line[0] == '#':
			continue
		p = breakfile(line)
		if 'hsa' not in p[1]:
			continue
		items = p[1]
		family2members.setdefault(p[0],[]).append(items)
		member2family_name[items] = p[0]








	return mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age

def disease_number_correlations(mirna2disease,mirna2age):
	age = []
	number_disease = []

	for mirna in mirna2disease:
		age.append(mirna2age[mirna])
		number_disease.append(len(mirna2disease[mirna]))

	final_corr = spearmanr(age,number_disease)
	fle = open('txtfles/corrs.txt','a')
	fle.write(str(final_corr) + '\n')
	fle.close()

	age2num = {}
	for mirna in mirna2disease:
		age2num.setdefault(mirna2age[mirna], []).append(len(mirna2disease[mirna]))

	labels = sorted(age2num.keys())
	nums = [age2num[i] for i in labels]
	plt.figure(figsize=(10,7))
	plt.boxplot(nums,labels=labels)
	plt.xlabel('miRNA Ages')
	plt.ylabel('Number of Associated Diseases')
	plt.title('miRNA Age versus Number of Disease Associations')
	plt.savefig('images/mirna_ages_vs_num_dis.png')
	plt.close()

def disease_txt_files(mirna2disease, disease2mirna):
	for disease in disease2mirna:
		name = disease.split(' ')
		name = '-'.join(name)
		newname = 'disease_mirnas/%s.txt' %(name)
		fle = open(newname,'w')

		for mirna in disease2mirna[disease]:
			fle.write('%s\n' %(mirna))
		fle.close()


def make_vector(mirna_name,mirna2disease, diseaselst):
	bin_vec = []

	for dis in diseaselst:
		if dis in mirna2disease[mirna_name]:
			bin_vec.append(str(1))
		else:
			bin_vec.append(str(0))
	return bin_vec



def hamming_distance(mirna2age, family2members, member2family_name,diseaselst, mirna2disease):
	mirna2bin_vec = {}
	family2disease_members = {}


	for fam in family2members:
		new_mems = []
		for mirna in family2members[fam]:
			if mirna in mirna2disease and mirna in mirna2age:
				new_mems.append(mirna)
		if len(new_mems) > 4:
			family2disease_members[fam] = new_mems


	mirna_lst = [alpha for x in family2disease_members.values() for alpha in x]

	for mirna in mirna_lst:
		mirna2bin_vec[mirna] = make_vector(mirna, mirna2disease, diseaselst)

	family2hamming_distances = {}
	mirna2family_hamming = {}


	for fam in family2disease_members:
		members = family2disease_members[fam]
		family_vec = []
		for mirna in members:
			mir_spec_vec = []
			for other in members:
				vec1 = ''.join(mirna2bin_vec[mirna])
				vec2 = ''.join(mirna2bin_vec[other])
				mir_spec_vec.append(hamming(vec1, vec2))
			mirna2family_hamming[mirna] = mir_spec_vec[:]
			family_vec.append(mir_spec_vec[:])
		family2hamming_distances[fam] = family_vec[:]

	family_average_age = []
	mirna_max_hamming = []

	for mirna in mirna2family_hamming:
		fam_name = member2family_name[mirna]	
		family_average_age.append(mean([mirna2age[mirna] for mirna in family2disease_members[fam_name]]))
		mirna_max_hamming.append( max(mirna2family_hamming[mirna]))

	mirna_age2hamming = {}

	for mirna in mirna2family_hamming:
		mirna_age2hamming.setdefault(float(mirna2age[mirna]),[]).append(max(mirna2family_hamming[mirna]))


	corr =  spearmanr(family_average_age, mirna_max_hamming)
	# print corr
	labels = sorted(mirna_age2hamming.keys())
	nums = [mirna_age2hamming[i] for i in labels]
	fig, ax1 = plt.subplots(figsize=(10,10))
	ax1.set_ylim(0, 75)
	plt.boxplot(nums,labels=labels)
	plt.xlabel('miRNA Family Age')
	plt.ylabel('Hamming Vector')
	plt.title('An Illustration of miRNA Age Versus the Similarity of Diseases It Targets')
	plt.savefig('images/mirna_ages_vs_hamming.png')
	plt.close()



def break_target(fle, tipo):
	fle = open(fle,'r')
	text = fle.readlines()
	fle.close()

	mirna2targets = {}
	targets2mirna = {}

	if tipo == 'verified':
		for line in text:
			p = breakfile(line)
			if '3p' in p[1] or '5p' in p[1]:
				mir = p[1][:-3].lower()
			else:
				mir = p[1].lower()
			target = p[3]
			if target == '':
				continue

			mirna2targets.setdefault(mir,[]).append(target)
			targets2mirna.setdefault(target,[]).append(mir)

		verified_dicts = [mirna2targets, targets2mirna]


	
	return verified_dicts




def target_mirna_corrs(verified_dicts,mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age):
	mirs_with_tar = []
	mirna2targets = verified_dicts[0]
	targets2mirna = verified_dicts[1]
	mirs_with_tar = mirna2targets.keys()

	mir_ages = []
	tar_ages = []
	for mirna in mirna2targets:
		for target in mirna2targets[mirna]:
			if target in gene2age and mirna in mirna2age:
				mir_ages.append(float(mirna2age[mirna]))
				tar_ages.append(float(gene2age[target]))




	mir_ages = []
	tar_nums = []

	for mirna in mirna2targets:
		if mirna in mirna2age:
			mir_ages.append(float(mirna2age[mirna]))
			tar_nums.append(len(mirna2targets[mirna]))
	# print spearmanr(mir_ages, tar_nums)

	age2num = {}
	for mirna in mirna2targets:
		if mirna not in mirna2age:
			continue
		age2num.setdefault(mirna2age[mirna], []).append(len(mirna2targets[mirna]))

	labels = sorted(age2num.keys())
	nums = [age2num[i] for i in labels]
	fig, ax1 = plt.subplots(figsize=(10,7))
	ax1.set_ylim(0, 50)
	plt.boxplot(nums,labels=labels)
	plt.xlabel('miRNA Ages')
	plt.ylabel('Number of Targets')
	plt.title('miRNA Age versus Number of Targets')
	plt.savefig('images/mirna_ages_vs_num_tars.png')
	plt.close()

	for mirna in mirna2targets:
		fle = open('txtfles/mirna2target/' + mirna + '.txt','w')
		for item in mirna2targets[mirna]:
			fle.write(item + '\n')
		fle.close()
	tarlst = []

	for mirna in mirna2targets:
		tarlst = tarlst[:] + mirna2targets[mirna]
	tarlst = list(set(tarlst))

	fle = open('txtfles/tarlst.txt','w')
	for item in tarlst:
		fle.write(item + '\n')
	fle.close()

	gene_corr_ages = []
	mir_ages = []
	for tar in targets2mirna:
		if tar in gene2age:
			ages = [mirna2age[i] for i in targets2mirna[tar] if i in mirna2age]
			if len(ages) == 0:
				continue
			for item in ages:
				mir_ages.append(float(item))
				gene_corr_ages.append(float(gene2age[tar]))

	boxplot(mir_ages, gene_corr_ages, 'miRNA Age (MYA) vs Gene Age (MYA)', 'miRNA Age', 'Gene Agae', 'mirna_age_vs_gene_age')

	mir_over_tar = []
	number_over = 0
	total = 0

	tarlst = list(set(targets2mirna.keys()))
	
	values4avg = []
	for mirna in mirna2targets:
		if mirna not in mirna2age:
			continue
		num = 0
		tars = [i for i in mirna2targets[mirna] if i in gene2age]
		if len(tars) == 0:
			continue
		mirage = float(mirna2age[mirna])
		for tar in tars:
			if float(gene2age[tar]) < mirage:
			 		num += 1
		values4avg.append(float(num)/ float(len(tars)))
	avg_under = mean(values4avg)
	# print avg_under


	interac_under = 0
	total_interac = 0

	for tar in targets2mirna:
		if tar in gene2age:
			mirna_ages = [float(mirna2age[i]) for i in targets2mirna[tar] if i in mirna2age]
			interac_under += int(len([i for i in mirna_ages if i >= float(gene2age[tar])]))
			total_interac += len(targets2mirna[tar])

	print float(interac_under) / float(total_interac)





	num = 0

	for i in range(1000):
		interac_under_ran = 0
		total_interac_ran = 0
		for mirna in mirna2targets:
			if mirna in mirna2age:
				tars = [i for i in mirna2targets[mirna] if i in gene2age]
				total_interac_ran += len(tars)
				ran_samp = [i for i in random.sample(gene2age.keys(),len(tars)) if float(gene2age[i]) <= mirna2age[mirna]]
				interac_under_ran += len(ran_samp)
		print float(interac_under_ran) / float(total_interac_ran) 
		if float(interac_under_ran) / float(total_interac_ran) <= .16:
			num += 1
	print num

	bincount('Target Ages', 'Ages of Targets', [float(gene2age[i]) for i in tarlst if i in gene2age], 'images/tar_ages.png')
	secondlst = []
	for tar in targets2mirna:
		if tar in gene2age:
			mirna_ages = [float(mirna2age[i]) for i in targets2mirna[tar] if i in mirna2age]
			if  len([i for i in mirna_ages if i >= float(gene2age[tar])]) != 0:
				secondlst.append(tar)
	bincount('Target Ages', 'Ages of Targets with Older miRNA Target ', [float(gene2age[i]) for i in secondlst], 'images/tar_ages_younger.png')



	percent_under_lst = []


	for mirna in mirna2targets:
		if mirna not in mirna2age:
			continue
		tar_ages = [float(gene2age[i]) for i in mirna2targets[mirna] if i in gene2age]
		if len(tar_ages) == 0:
			continue
		percent = 0.0
		counter_under = 0
		for i in tar_ages:
			if i <= float(mirna2age[mirna]):
				counter_under += 1
		percent_under_lst.append(float(counter_under) / float(len(tar_ages)))

	a = np.array(percent_under_lst[:])

	plt.hist(a,20)
	plt.xlabel('Percent of Targets Under The miRNA Age')
	plt.ylabel('Frequency')
	plt.title('Distribution of Percent of Targets Under The miRNA Age')
	plt.savefig('images/percent_under_hist.png')
	
	plt.close()












def enrichment_lists(verified_dicts,mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age):
	fle = open('mirnas_in_disease.txt','w')

	for mirna in mirna2age:
		if mirna in mirna2disease:
			fle.write(mirna + '\n')

	fle.close()

	fle = open('mirnas_not_in_disease.txt','w')

	for mirna in mirna2age:
		if mirna not in mirna2disease:
			fle.write(mirna + '\n')
	fle.close()


	diseases = disease2mirna.keys()
	fle = open('txtfles/disease_lst.txt','w')

	for dis in diseases:
		fle.write(dis + '\n')

	fle.close()

	disease2stats = {}


	for dis in diseases:

		fle = open('/Users/virpatel/projects/vanderbilt-summer-2014/main_script/disease_mirnas_ages/mirbase_method2_ages--' + string.replace(dis,' ', '-')+ '--Background_summary.txt')
		text = fle.readlines()
		mirnalst = []

		median = 0
		avg = 0
		mann_value_andp = []

		for line in text:
			
			if line[0] == '#':
				if 'Mann-Whitney U test: U =' not in line and 'Average age of' not in line and ' Median age of' not in line :
					continue
			if line[0] != '#':
				p = line.split('\t')
				mirnalst.append([float(p[0]),int(p[2]),int(p[3]),float(string.replace(p[4],'\n',''))])
			else:
				if 'Median' in line:
					if 'Background' in line:
						continue
					p = line.split(':')
					median = float(string.replace(p[1],' ',''))
				if 'Average' in line:
					if 'Background' in line:
						continue
					p = line.split(':')
					avg = float(string.replace(p[1],' ', ''))

				if 'Mann-Whitney U test: U' in line:

					part1 = string.replace(line.split('U = ')[1],')\n','')
					part2 = part1.split(' (p = ')
					mann_value_andp = part2[:]


					


		biglst = [dis,mann_value_andp, median, avg, mirnalst]
		disease2stats[dis] = biglst


		fle.close()

	fle = open('txtfles/diseases_enrich.txt','w')
	fle.write('Name\tAverage Age\tMann-Whitney U Value\tNumber of miRNAs\n')
	for item in disease2stats:
		a = disease2stats[item]
		lst = a[-1]
		total = 0
		for el in lst:
			total += el[1]
		fle.write('%s\t%.1f\t%s\t%i\n' %(a[0], a[3], a[1][0], total))
	fle.close()

	biggest = ['Inflammation']
	for dis in diseases:
		if disease2stats[biggest[0]][3] > disease2stats[dis][3]:
			continue
		elif disease2stats[biggest[0]][3] < disease2stats[dis][3]:
			biggest = [dis]
		elif disease2stats[biggest[0]][3] == disease2stats[dis][3]:
			biggest.append(dis)


	smallest = ['Inflammation']
	for dis in diseases:
		if disease2stats[smallest[0]][3] < disease2stats[dis][3]:
			continue
		elif disease2stats[smallest[0]][3] > disease2stats[dis][3]:
			smallest = [dis]
		elif disease2stats[smallest[0]][3] == disease2stats[dis][3]:
			smallest.append(dis)

	fle = open('txtfles/disease2avg.txt','w')
	disease2avgage = {}
	for dis in diseases:
		age = disease2stats[dis][3]
		disease2avgage.setdefault(age,[]).append(dis)

	
	for age in sorted(disease2avgage):
		fle.write('%s\t%s\n' % (age, str(disease2avgage[age])))
	fle.close()



def stability_test(verified_dicts, mirna2age, age2mirna, disease2mirna, mirna2disease, age2disease, disease2age, family2members, member2family_name, gene2age, stab_text):
	fle = open(stab_text, 'r')
	txt = fle.readlines()
	fle.close()

	line2mirnas = {}
	mirnas2stab = {}
	cell_lines = breakfile(txt[0])[1:]
	for index,cell in enumerate(cell_lines):
		mirnas2stab = {}
		for line in txt[1:]:
			if line[0:3] != 'hsa':
				continue
			else:
				p = breakfile(line)
				if p[index + 1] == 'N/A':
					continue
				if '3p' in p[0] or '5p' in p[0]:
					name = p[0].split('-')
					name = '-'.join(name[:-1]).lower()
					mirnas2stab[name] = float(p[index + 1])
					continue
				mirnas2stab[p[0].lower()] = float(p[index + 1])
		line2mirnas[cell] = mirnas2stab
	
	for cell in line2mirnas:
		ages = []
		stab = []
		counter = 0
		for mirna in line2mirnas[cell]:
			if mirna not in mirna2age:
				continue
			counter += 1
			ages.append(float(mirna2age[mirna]))
			stab.append(float(line2mirnas[cell][mirna]))
		bins_4_pic = {}
		labels = sorted(list(set(ages)))
		for index, num in enumerate(stab):
			bins_4_pic.setdefault(ages[index],[]).append(num)
	nums = [bins_4_pic[i] for i in labels]
	plt.figure(figsize=(10,7))
	plt.boxplot(nums,labels=labels)
	plt.xlabel('miRNA Ages')
	plt.ylabel('miRNA Stability')
	plt.title('miRNA Age versus MSI')
	plt.savefig('images/mirna_ages_vs_stability1.png')
	plt.close()





def utr_stuff(verified_dicts, mirna2age, age2mirna, disease2mirna, mirna2disease, age2disease, disease2age, family2members, member2family_name, gene2age,utr_txt):
	fle = open('txtfles/3utr.txt','r')
	txt = fle.readlines()
	fle.close()

	gene2utr = {}

	for line in txt:
		p = breakfile(line)
		length = abs(int(p[1]) - int(p[2]))
		new =  p[3].split('_')
		new = new[0] + '_' + new[1]
		gene2utr.setdefault(new,[]).append(length)


	fle = open('txtfles/hgnc_stuff.txt','r')
	txt = fle.readlines()
	fle.close()

	nm2hgnc = {}

	for line in txt:
		p = breakfile(line)
		nm2hgnc[p[-1]] = p[1]


	hgnc2utr = {}

	for gene in gene2utr:
		if gene in nm2hgnc:
			hgnc2utr[nm2hgnc[gene]] = gene2utr[gene]

	gene2utr = {}
	fle = open('txtfles/gene2avg_utr.txt','w')
	for gene in hgnc2utr:
		length = np.mean(hgnc2utr[gene])
		gene2utr[gene] = float(length)

	fle.close()

	mirna2targets = verified_dicts[0]
	targets2mirna = verified_dicts[1]
	mirna2utr_lengths = {}

	for mirna in mirna2targets:
		for target in mirna2targets[mirna]:
			if target in gene2utr:
				mirna2utr_lengths.setdefault(mirna,[]).append(gene2utr[target])

	ages = []
	utr_lens = []

	for mirna in mirna2utr_lengths:
		if mirna in mirna2age:
			# print mirna2utr_lengths[mirna]
			avg = np.mean(mirna2utr_lengths[mirna])
			utr_lens.append(avg)
			ages.append(float(mirna2age[mirna]))

	labels = sorted(list(set(ages)))
	bin4vec = {}
	for index, item in enumerate(utr_lens):
		bin4vec.setdefault(ages[index],[]).append(item)
	nums = [bin4vec[i] for i in labels]
	plt.figure(figsize=(10,7))
	plt.boxplot(nums,labels=labels)
	plt.xlabel('miRNA Ages')
	plt.ylabel('UTR Lengths')
	plt.title('miRNA Age versus UTR Lengths')
	plt.savefig('images/mirna_ages_vs_utr.png')
	plt.close()


	











	





def main():


	mirna_ages = ''
	gene_ages = ''
	disease_associations = ''
	family_associations = ''
	predicted_targets = ''
	verified_targets = ''



	args = sys.argv

	for index, arg in enumerate(args):
		if arg == '-age_file':
			mirna_ages = args[index + 1]
			continue
		if arg == '-gene_ages':
			gene_ages = args[index + 1]
			continue
		if arg == '-diseases':
			disease_associations = args[index + 1]
			continue
		if arg == '-families':
			family_associations = args[index + 1]
		if arg == '-verified_targets':
			verified_targets = args[index + 1]
		if arg == '-predicted_targets':
			predicted_targets = args[index + 1]


	fle = open('txtfles/corrs.txt',"w")
	fle.close()


			

	mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age = break_files(mirna_ages, gene_ages, disease_associations, family_associations)

	disease_number_correlations(mirna2disease, mirna2age)




	a = [dis for alpha in mirna2disease.values() for dis in alpha]
	diseaselst = sorted(list(set(a)))

	family2hamming_distances = hamming_distance(mirna2age, family2members, member2family_name, diseaselst, mirna2disease)
	mirna_rates(mirna2age)
	verified_dicts = break_target(verified_targets, 'verified')

	target_mirna_corrs(verified_dicts,mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age)

	enrichment_lists(verified_dicts, mirna2age, age2mirna, disease2mirna, mirna2disease, age2disease, disease2age, family2members, member2family_name, gene2age)
	disease_txt_files(mirna2disease, disease2mirna)

	utr_stuff(verified_dicts, mirna2age, age2mirna, disease2mirna, mirna2disease, age2disease, disease2age, family2members, member2family_name, gene2age, 'txtfles/3utr.txt')

	stability_test(verified_dicts, mirna2age, age2mirna, disease2mirna, mirna2disease, age2disease, disease2age, family2members, member2family_name, gene2age, 'txtfles/MSI.txt')

main()



