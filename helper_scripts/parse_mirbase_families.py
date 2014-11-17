import filebreak
fle = open('mirbase_families.txt','r')
text = fle.readlines()
fle.close()

family_2_members = {}

current_fam = "vir"
for line in text:
	if line[:2] == "ID":
		pie = filebreak.breakfile(line)
		current_fam = pie[1].lower()
		continue
	if line[:2] == "AC" or line[:2] == "//":
		continue
	else:
		pie = filebreak.breakfile(line)
		family_2_members.setdefault(current_fam,[]).append(pie[2].lower())

fle = open('mirbase_parsed_families.txt',"w")

for fam in family_2_members:
	for item in family_2_members[fam]:
		fle.write('%s\t%s\n' % (fam,item))

fle.close()