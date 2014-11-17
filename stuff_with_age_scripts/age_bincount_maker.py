# general file that will create a bincount for any age file from PH

import filebreak
import virplot
fle = open(ANY FILE YOU WISH TO COUNT,'r')
text = fle.readlines()
fle.close()
mirna2age = {}

for line in text:
	if line[0] == "#":
		continue
	else:
		pie = filebreak.breakfile(line)
		mirna2age[pie[0]] = float(pie[1])

age_list_master = mirna2age.values()
virplot.bincount('miRNA Ages', 'Method 1 mirBase Ages', age_list_master, 'method1_mirbase_bin.png') # replace argument 1 with name of graph, argument 2 with x axis label, argument 4 with name of image along with file type

