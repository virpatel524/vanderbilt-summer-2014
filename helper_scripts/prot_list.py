import filebreak


fle = open("prot.txt")


text = fle.readlines()
fle.close()

fle = open("hgnc.txt","w")

for index,line in enumerate(text):
	if index == 0:
		continue
	else:
		prot = filebreak.breakfile(line)
		prot = prot[1]
		fle.write(prot + "\n")


fle.close()