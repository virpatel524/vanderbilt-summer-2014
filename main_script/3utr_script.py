import string
import re
import scipy
import matplotlib.pyplot as plt
import calendar
import time
import sys, os
import math
import random
import numpy as np


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


fle = open('txtfles/gene2avg_utr.txt','w')
for gene in hgnc2utr:
	length = np.mean(hgnc2utr[gene])
	fle.write('%s\t%d\n' %(gene,int(length)))

fle.close()