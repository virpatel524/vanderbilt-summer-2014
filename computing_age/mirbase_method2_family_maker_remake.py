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





fle = open("organisms.txt","r")
text = fle.readlines()
fle.close()

code2rest = {}
for line in text:
	if line[0] == '#':
		continue
	p = breakfile(line)
	code2rest[p[0]] = p[3]


bad_ones = ['aqc','ghb','ssp','fru','hsv','hvt','ebv','rlc','hhv','mcm','pbi','jcv','bkv','mdv','hma','bpc','ksh']



fle = open("mirbase_parsed_families.txt","r")
text = fle.readlines()
fle.close()

fam2members = {}

for line in text:
	p = breakfile(line)
	code = p[1].split('-')[0]
	if 'Virus' in code2rest[code]:
		continue
	if code in bad_ones:
		continue
	fam2members.setdefault(p[0],[]).append(p[1])

fle = open('famfilenew.txt','w')

for fam in fam2members:
	stringer = ''
	mirnas = fam2members[fam]
	for index, mirna in enumerate(mirnas):
		code = mirna.split('-')[0]
		if mirna != mirnas[-1]:
			stringer += "%s|mirBase:%s " % (code,mirna)
		else:
			stringer += "%s|mirBase:%s" % (code,mirna)
	fle.write(stringer + "\n")


fle.close()
