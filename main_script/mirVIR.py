import sys, os, string, numpy, matplotlib.pyplot as plt, basicmethods

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





def break_files(mirna_ages,gene_ages,diseases,tree):
	fle1 = open(mirna_ages,'r')
	mirna_ages_unparsed = fle1.readlines()
	fle1.close()

	fle2 = open(gene_ages,'r')
	gene_ages_unparsed = fle2.readlines()
	fle2.close()

	fle3 = open(diseases,'r')
	disease_associations = fle3.readlines()
	fle3.close()

	fle4 = open(tree,'r')
	tree_text = fle4.readlines()[0]
	fle4.close()



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


	return mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age

def main():


	mirna_ages = ''
	gene_ages = ''
	disease_associations = ''
	family_associations = ''


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
		if arg == '-tree' in arg:
			family_associations = args[index + 1]
			

	mirna2age,age2mirna,disease2mirna,mirna2disease,age2disease, disease2age = break_files(mirna_ages, gene_ages, disease_associations, family_associations)

	




main()



