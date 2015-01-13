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


def bincountdict(lst):
	countdict = {}
	for item in lst:
		if item in countdict:
			countdict[item] += 1
		else:
			countdict[item] = 1

	return countdict


fle = open("mirbase_original_data.txt","r")
text = fle.readlines()
fle.close()
mirna2mem = {}

for line in text:
	if line[:2] != "MI":
		continue
	else:
		pie = breakfile(line)
		print pie
		part1 = pie[2][4:]
		part2 = pie[2][:3]
		mirna2mem.setdefault(part1,[]).append(part2)


fle = open("organisms.txt","r")
text = fle.readlines()
fle.close()
ncbi2name = {}
species2name = {}
name2both = {}
for line in text:
	pie = breakfile(line)
	if "Virus" in pie[3]:
		continue
	ncbi2name[pie[4]] = pie[0]
	species2name[pie[2]] = pie[0]
	name2both[pie[0]] = [pie[2],pie[4]]

poplst = ["ssp","ghb"]
for el in poplst:
	name2both.pop(el,None)

badapples = []

for mirna in mirna2mem:
	mirna2mem[mirna] = [mem for mem in mirna2mem[mirna] if mem in name2both.keys()]
	if mirna2mem[mirna] == []:
		badapples.append(mirna)

for item in badapples:
	mirna2mem.pop(item,None)


fle = open("famfilenew.txt","w")
bad_ones = ['aqc','ghb','ssp','fru','hsv','hvt','ebv','rlc','hhv','mcm','pbi','jcv','bkv','mdv','hma','bpc','ksh']
for mirna_upper in mirna2mem:
	mirna = mirna_upper.lower()
	species = [mir for mir in mirna2mem[mirna_upper] if mir not in bad_ones]
	if species == []:
		continue
	string = ""
	for index, item in enumerate(species):
		if item != species[-1]:
			string += "%s|mirBase:%s-%s " % (item,item,mirna)
		else:
			string += "%s|mirBase:%s-%s" % (item,item,mirna)
	fle.write(string + "\n")



fle.close()




