# easily extendable script that does all of my analyses
import numpy as np
from numpy.random import randn
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import string
import re
import scipy
import calendar
import time
import sys, os
import math
import random


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




def list_write(fle,lst,para="t"):
	for index,el in enumerate(lst):
		if para == "t":
			if index != len(lst) - 1:
				fle.write(str(el) + "\t")
			else:
				fle.write(str(el) + "\n")
		if para == "n":
			fle.write(str(el) + "\n")


def bincount(xlab,name,lst,filename):
	data = lst
	alldata = sorted(data)
	labs = sorted(list(set(data)))
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


def mirna_rates(mirna2age):
	mirna_ages = sorted(list(set(mirna2age.values())))
	number_in_agecount = {}
	for age in mirna_ages:
		number_in_agecount[age] = 0

	for age in mirna2age.values():
		number_in_agecount[age] = number_in_agecount[age] + 1

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
		gene2age[p[0]] = p[1]







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
				diseaselst = diseaselst + mirna2disease[mirna]
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
	print final_corr

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
	plt.savefig('mirna_ages_vs_num_dis.png')
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
			mirna2family_hamming[mirna] = mir_spec_vec
			family_vec.append(mir_spec_vec)
		family2hamming_distances[fam] = family_vec

	family_average_age = []
	mirna_max_hamming = []

	for mirna in mirna2family_hamming:
		fam_name = member2family_name[mirna]	
		family_average_age.append(mean([mirna2age[mirna] for mirna in family2disease_members[fam_name]]))
		mirna_max_hamming.append( max(mirna2family_hamming[mirna]))

	corr =  spearmanr(family_average_age, mirna_max_hamming)
	print corr



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

			mirna2targets.setdefault(mir,[]).append(target)
			targets2mirna.setdefault(target,[]).append(mir)

		verified_dicts = [mirna2targets, targets2mirna]


	
	return verified_dicts




def target_mirna_corrs(verified_dicts,mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age):
	mirs_with_tar = []
	mirna2targets = verified_dicts[0]
	targets2mirna = verified_dicts[1]
	for mirna in mirna2age:
		if mirna in mirna2targets:
			mirs_with_tar.append(mirna)

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
	print spearmanr(mir_ages, tar_nums)

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
	plt.savefig('mirna_ages_vs_num_tars.png')
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


			

	mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age = break_files(mirna_ages, gene_ages, disease_associations, family_associations)


	disease_number_correlations(mirna2disease, mirna2age)




	a = [dis for alpha in mirna2disease.values() for dis in alpha]
	diseaselst = sorted(list(set(a)))
	family2hamming_distances = hamming_distance(mirna2age, family2members, member2family_name, diseaselst, mirna2disease)
	mirna_rates(mirna2age)
	verified_dicts = break_target(verified_targets, 'verified')

	target_mirna_corrs(verified_dicts,mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age, family2members, member2family_name, gene2age)


	disease_txt_files(mirna2disease, disease2mirna)



main()



