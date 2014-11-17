import scripter

fle = open("mirname2species.txt")
text1 = fle.readlines()
fle.close()

fle = open("bigtreewitplants.txt")
text2 = fle.readlines()
fle.close()

mir2species = {}

for line in text1:
	pie = scripter.basicmethods.breakfile(line)
	mir2species[pie[0]] = pie[1]

tree = text2[0]

missing = []

for item in mir2species:
		treecopy = tree[:]
		old = scripter.string.replace(mir2species[item], " ", "_")
		tree = scripter.string.replace(tree, old, item)
		if tree == treecopy:
			missing.append(item)

fle = open("newtree.txt","w")
fle.write(tree)
fle.close()