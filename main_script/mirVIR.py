import sys, os, string, numpy, matplotlib.pyplot as plt, basicmethods

def breakfile(line):
	temp = line
	temp = string.replace(temp, "\n", "")
	temp = temp.split("\t")

	return tem


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
	tree_file = fle4.readlines()
	fle4.close()








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
		

break_files(mirna_ages, gene_ages, disease_associations, family_associations)






