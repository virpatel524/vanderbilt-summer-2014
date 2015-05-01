import basicmethods, string

fle = open('final_tree_raw.txt','r')
oldtree = fle.readlines()[0]
fle.close()


fle = open('three_letter_species.txt','r')
species_list = fle.readlines()
fle.close()

letter2name = {}
name2letter = {}

for line in species_list:
	p = basicmethods.breakfile(line)
	letter2name[p[0]] = p[1]
	name2letter[p[1]] = p[0]

for key in letter2name:
	replacer = string.replace(letter2name[key],' ', '_')
	if replacer in oldtree:
		oldtree = string.replace(oldtree,replacer,key)

fle = open('final_tree_replaced.txt','w')
fle.write(oldtree)
fle.close()




