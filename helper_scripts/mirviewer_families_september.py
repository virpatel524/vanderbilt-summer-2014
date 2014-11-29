# parser of mirna families

import basicmethods

fle = open(# file we want to put in,"r")
text = fle.readlines()
fle.close()

mirnalst = []

for line in text:
	pie = basicmethods.breakfile(line)[0]
	pie  = pie.split(" ")
	for el in pie:
		allstuff = el.split(":")
		mirnalst.append(el[1])
		

fle = open("family_2_members_mirviewer.txt","r")
text = fle.readlines()
fle.close()

familydict = {}

for line in text:
	if line[0] == "#":
		continue
	