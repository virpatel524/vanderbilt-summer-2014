import matplotlib.pyplot as plt,string


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




fle = open('mirbase_method2_ages.txt','r')
ages = fle.readlines()
fle.close()
ageslst = []
for line in ages:
	if line[0] == '#':
		continue
	p = breakfile(line)
	ageslst.append(float(p[1]))

bincount('Ages', 'Bincount of Ages in mirBase Method 2' ,ageslst, 'mirbase_method2_ages_bincount.png')


	



