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
import matplotlib.pyplot as plt
import calendar
import time
import sys, os
import math
import random



def breakfile(line):
	temp = line
	temp = string.replace(temp, "\n", "")
	temp = temp.split("\t")

	return temp

def bincount(xlab,name,lst,filename):
	data = [float(a) for a in lst]
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



arrg = sys.argv
gene_fle = ''
for index, arg in enumerate(arrg):
	if arg == '-age_file':
		gene_fle = arrg[index + 1]


fle = open(gene_fle, 'r')
text = fle.readlines()
fle.close()

gene2age = {}

for line in text:
	p = breakfile(line)
	gene2age[p[0]] = p[1]


age_lst = gene2age.values()


bincount('Ages', 'Bincount of Gene Ages (PH)' ,age_lst, 'ph_ages_bincount.png')







