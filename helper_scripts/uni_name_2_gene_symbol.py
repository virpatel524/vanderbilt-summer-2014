# a helper script. not that important. used to change the names of the proteins 

import string

def breakfile(line):
	temp = line
	temp = string.replace(temp, "\n", "")
	temp = temp.split("\t")

	return temp






fle1 = open("wagner.txt","r")
text = fle1.readlines()

fle2 = open("dollo.txt","r")
text2 = fle2.readlines()

uni2age = {}
uni2sym = {}
sym2age = {}

for line in text:
	if line[0] == "#":
		continue

	pie = breakfile(line)
	try:
		var = float(pie[2])
	except ValueError:
		var = pie[2]
	print pie[2]

	if isinstance(var, float):
		continue

	else:
		uni2sym[pie[1]] = pie[2]

for line in text2:
	if line[0] == "#":
		continue

	pie = breakfile(line)

	uni2age[pie[1]] = pie[2]


print uni2sym

for uni in uni2age:
	if uni in uni2sym:
		sym2age[ uni2sym[uni]] = uni2age[uni]


finfle = open("sym2age.txt","w")

for sym in sym2age:
	finfle.write("%s\t%s\n" %(sym, sym2age[sym]))













