
# generates family files for mirbase data using method 1

import filebreak
fle = open('mirbase_families.txt','r')
text = fle.readlines()
fle.close()
mirna2species = {}
for line in text:
	if line[:2] != "MI":
		continue
	pie = filebreak.breakfile(line)
	full_mirna = pie[2].lower()
	species = full_mirna[:3]
	if species in ['ghb','ssp','fru','hsv','hvt','ebv','rlc','hhv','mcm','pbi','jcv','bkv','mdv','hma','bpc','ksh']:
		continue
	if species[-1] == 'v':
		continue
	mirna_name = full_mirna[4:]
	mirna2species.setdefault(mirna_name,[]).append(species)

fle = open('mirbase_method1_def.txt','w')

for mir in mirna2species:
	s = ''
	for index,species in enumerate(mirna2species[mir]):
		if index == len(mirna2species[mir]) - 1:
			s += '%s|mirBase:%s-%s\n' %(species,species,mir)
		else:
			s += '%s|mirBase:%s-%s ' %(species,species,mir)

	fle.write(s)