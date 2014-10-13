import string
import time
import matplotlib.pyplot as plt


def getdate():
	return (time.strftime("%m/%d/%Y"))


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


