# computes rates in new miRNAs / MYA for each clade divison per TimeTree estimates
import filebreak, sys
agefle = ''
for index, arg in enumerate(sys.argv):
	if arg = '-agefile': agefle = sys.argv[index + 1] 

fle = open(agefle,"r")
text = fle.readlines()
ages = {}
for line in text:
	if line[0] == "#":
		continue
	else:
		pie = filebreak.breakfile(line)
		temp = float(pie[1])

		if temp in ages:
			ages[temp] += 1
		else:
			ages[temp] = 1


keys = sorted(ages.keys())

relative = []

for index,item in enumerate(keys):
	if item == keys[-1]:
		continue
	per = float(ages[item]) / float(keys[index + 1] - keys[index])
	relative.append(per)

print relative

print 72/2 *100











