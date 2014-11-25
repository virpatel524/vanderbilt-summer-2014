import scripter

fle = open("miFam.dat","r")
text = fle.readlines()
fle.close()
mirna2mem = {}

for line in text:
	if line[:2] != "MI":
		continue
	else:
		pie = scripter.basicmethods.breakfile(line)
		print pie
		part1 = pie[2][4:]
		part2 = pie[2][:3]
		mirna2mem.setdefault(part1,[]).append(part2)


fle = open("organ.txt","r")
text = fle.readlines()
fle.close()
ncbi2name = {}
species2name = {}
name2both = {}
for line in text:
	pie = scripter.basicmethods.breakfile(line)
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
for mirna in mirna2mem:
	species = mirna2mem[mirna]
	string = ""
	for index, item in enumerate(species):
		if item == "ssp" or item == "ghb":
			continue
		if item != species[-1]:
			string += "%s|mirBase:%s-%s " % (item,item,mirna)
		else:
			string += "%s|mirBase:%s-%s" % (item,item,mirna)
	fle.write(string + "\n")



fle.close()




