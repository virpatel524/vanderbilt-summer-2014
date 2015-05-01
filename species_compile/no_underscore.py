import basicmethods,string

# mirbase

fle = open('mirbase_abbreviations.txt','r')
txt = fle.readlines()
fle.close()

fle = open('mirbase_species.txt','w')

for line in txt:
	p = basicmethods.breakfile(line)
	three_letter = p[0]
	ncbi = p[4]
	species = p[2]
	if 'Virus' in p[3]:
		continue
	if '#' in line:
		continue
	fle.write('%s\t%s\t%s\n' % (three_letter, species, ncbi))

fle.close()


# mirviewer

fle = open('species_underscore.txt','r')
txt = fle.readlines()
fle.close()

fle = open('mirviewer_species.txt','w')

for line in txt:
	p = basicmethods.breakfile(line)
	three_letter = p[0]
	species = string.replace(p[1], '_', ' ')
	fle.write('%s\t%s\n' % (three_letter,species ))

fle.close()


# combine

fle = open('mirviewer_species.txt','r')
mirviewer = fle.readlines()
fle.close()

fle = open('mirbase_species.txt','r')
mirbase = fle.readlines()
fle.close()

species_list = []


for line in mirviewer + mirbase:
	p = basicmethods.breakfile(line)
	if p[1] not in species_list:
		species_list.append(p[1])



fle = open('combined_species_list_phylot.txt','w')
for item in species_list:
	if item not in ['Xenopus tropicalis']:
		fle.write(item + '\n')
	else:
		if 'Xenopus' in item:
			fle.write('8364\n')


fle.close()


fle = open('three_letter_species.txt','w')

letter_item_dict = {}

for item in mirviewer + mirbase:
	p = basicmethods.breakfile(item)
	if p[1] not in letter_item_dict.values():
		if p[0] in letter_item_dict.keys(): print p[0]
		letter_item_dict[p[0]] = p[1]

for key in letter_item_dict:
	fle.write('%s\t%s\n' %(key, letter_item_dict[key]))
fle.close()
