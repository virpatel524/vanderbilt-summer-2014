# long script used to compute most of my inital analyses.
# much of it is commented out 
__author__ = 'virpatel'


import itertools,string,numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import random
import math
import distance
import scipy.stats as scs
import scipy.cluster as cls
import operator
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

## vecompare method

def veccompare(vec1,vec2):
    total = len(vec1)
    sims = 0
    for index in range(0,len(vec1)):
        if vec1[index] == vec2[index]:
            sims+=1

    return distance.hamming(vec1, vec2, normalized=True)


# Dictionary Creation
################################################################################


fle = open("homomirna.txt","r")
holder = fle.readlines()
mirna=[]
datas= []
disease2mirna = {}
mirna2disease = {}
data = []
mirnaviewer = []
diseases = []
mirnaraw = []
nums = [] 

for i in holder:
    mirnaviewer.append(i.split("\n")[0])




fle = open("alldata.txt","r")
holder = []
holder = fle.readlines()
for i in holder:
    datas.append(i.split("\n")[0])



diseasedb = []
for var in datas:
    temp = var.split("\t")
    diseasedb.append([temp[1],temp[2]])

fle = open("hsa/hsa_95fams_dollo_age-depth.protein_list","r")
text = fle.readlines()
agesdb = []
for i in text:
    if i[0]!="#":
        temp = string.replace(i, "\n", "")
        temp = temp.split("\t")
        agesdb.append(temp)

agesdbages = []
agesdbmirna  = []
diseasedbdis = []
diseasedbmirna = []
mirnas = []


for i in diseasedb:
    diseasedbmirna.append(i[0])
    diseasedbdis.append(i[1])
for i in agesdb:
    agesdbmirna.append(i[0])
    agesdbages.append(i[1])

mirnacom = list(set(diseasedbmirna).intersection(set(agesdbmirna)))

mirna2age = {}
for i in agesdb:
    mirna2age.setdefault(i[0],i[1])


disease2mirna = {}
for i in diseasedb:
    if i[0] not in mirnacom: continue

    else:
        disease2mirna.setdefault(i[1],[]).append(i[0])


for key in disease2mirna:
    disease2mirna[key] = list(set(disease2mirna[key]))



disease2age = {}


for key in disease2mirna:
    newlst = []
    for i in disease2mirna[key]:
        newlst.append(mirna2age[i])
    disease2age[key] = newlst



badapples = []
for key in disease2age:
    if len(disease2age[key]) < 4:
        badapples.append(key)



for i in badapples:
    disease2age.pop(i,None)


for i in badapples:
    disease2mirna.pop(i,None)


age2mirna = {}
for i in agesdb:
    age2mirna.setdefault(i[1],[]).append(i[0])
fle.close()


fle = open("mirnasindisease.txt","w")
for i in mirnacom:

    fle.write(i+"\n")

fle.close()
fle = open("mirnanoindisease.txt","w")
for i in agesdb:
    contbool = True
    for key in disease2mirna:
        if i[0] in disease2mirna[key]: contbool=False
    if contbool == True:
        fle.write(i[0]+"\n")

# Bincounts

################################################################################

# ranger = range(0,17)
# for key in disease2age:
#
#     temp = key
#     temp = string.replace(key," ","")
#     stringer = ""
#     fle = open("diseaseMiRNAS/"+temp+"_mirnas.txt","w")
#     for i in disease2mirna[key]:
#         stringer+=i+"\n"
#     fle.write(stringer)
#     binlst = []
#     ints = disease2age[key]
#     for i in ints:
#         binlst.append(int(i))
#
#     binlst.sort()
#
#     freq = np.bincount(np.array(binlst),weights=None,minlength=17)
#
#
#
#
#
#     pos = np.arange(len(ranger))
#
#     width = 1.0
#     ax = plt.axes()
#
#     ax.set_xticks(pos + (width / 2))
#     ax.set_xticklabels(ranger)
#
#     plt.bar(pos,freq,width)
#     plt.xlim(pos.min(),pos.max()+width)
#
#     plt.savefig("binplotdisease/"+str(temp)+"_binplot.png")
#
#     plt.close()
#
#
#


# percentagedict = {}
#
#
# for key in disease2age:
#     stringer = disease2age[key]
#     ages = []
#     for i in stringer:
#         ages.append(int(i))
#
#     bins = [0]*17
#     total = len(ages)
#     percentages  = [0]*17
#     for i in ages:
#         bins[i]+=1
#     counter = 0
#     for i in range(0,len(percentages)):
#         percentages[i] = float(bins[i])/float(total)
#     percentagedict[key] = percentages
#


################################################################################

fle.close()


# Saturation Check

#
# for j in [10,15,20]:
#     for i in range(0,17):
#         fle = open("saturated"+str(j)+"/"+str(i)+"_"+str(j)+"percentsat.txt","w")
#
#         for key in percentagedict:
#             if percentagedict[key][i] >= float(j)/100.0:
#                 fle.write(key+"\n")
#         fle.close()
# #
# fle = open("cancerlist.txt","w")
# for key in disease2mirna:
#     fle.write(key+"\n")

#
# fle.close()
# fle = open("cancerlist.txt","r")
#
# cancers = fle.readlines()
# cancer = []
# for i in cancers:
#     cancer.append(string.replace(i,"\n",""))
# fle.close()
#
# cancermirna = []
#
# for i in cancer:
#     temp = disease2mirna[i]
#     for el in temp:
#         if el not in cancermirna and el in mirna2age:
#             cancermirna.append(el)
#


#
#
# fle = open("cancermirna.txt","w")
#
# for i in cancermirna:
#     fle.write(i+"\n")
#


# Dictionary Creation

################################################################################

for theta in mirna2age.keys():
    mirna2disease[theta] = []

for key in mirna2disease:
    for element in diseasedb:
        if element[0] == key and element[1] not in mirna2disease[key]:
            mirna2disease[key].append(element[1])



ageslst = []
numdiseases = []
for key in mirna2age:
    if mirna2disease[key] == [] : continue
    ageslst.append(mirna2age[key])
    numdiseases.append(len(mirna2disease[key]))

diseases = disease2mirna.keys()

disease2index = {}
for i in range(0,len(diseases)):
    disease2index[diseases[i]] = i

################################################################################

# Print a list of diseases

fle.close()
fle = open("textfiles/diseaselst.txt","w")

for d in diseases:
    fle.write(d+"\n")
fle.close()

mirna2vector = {}

# Create the vectors
################################################################################

def binvect(diseaselst):
    vect = [0]*len(diseases)
    for dis in diseaselst:
        if dis in disease2index:
            vect[disease2index[dis]] = 1

    return vect


def spebinvect(diseaselst):
    vect = [0]*len(spedis)
    for dis in diseaselst:
        if dis in spedisease2index :
            vect[spedisease2index[dis]] = 1

    return vect

################################################################################


for key in mirna2disease:
    mirna2vector[key] = (binvect(mirna2disease[key]))

file = open("fam2homo.txt","r")


text = file.readlines()


for index,oldstr in enumerate(text):
    text[index] = string.replace(oldstr,"\n","")
    text[index] = text[index].split("\t")
    text[index][0] = "hsa-"+text[index][0]


fam2kids = {}


for fam in text:
    fam2kids[fam[0]] = fam[1:]

for key in fam2kids:
    thekids = fam2kids[key]
    finalkids = []
    diseasemirna = mirna2disease.keys()
    for i in thekids:
        finalkids.append(i)
    fam2kids[key] = finalkids









badapples = []
for key in fam2kids:
    if fam2kids[key] == []: badapples.append(key)


for el in badapples:
    fam2kids.pop(el)

fam2kidvec = {}
for fam in fam2kids:
    davecs = []
    for el in fam2kids[fam]:
        if el in mirna2disease:
            davecs.append(binvect(mirna2disease[el]))
    fam2kidvec[fam] = davecs

similarities = {}




def ranvec():
    vec = []
    for i in range(0,len(diseases)):
        vec.append(random.randint(0,1))
    return vec
def vectable(table):
    biglst = []
    for i in table:
        minilst = []
        for j in table:
            minilst.append(veccompare(i,j))
        biglst.append(minilst)
    return biglst



similarities = {}


for key in fam2kidvec:
    similarities[key] = vectable(fam2kidvec[key])

fle.close()
for key in similarities:
    fle = open("sims4disease/"+key+".txt","w")

    for i in similarities[key]:
        fle.write(str(i) + "\n")
    fle.close()


def avgage(mirlst):
    agelst = []

    for mir in mirlst:
        if mir in agesdbmirna:
            agelst.append(float(mirna2age[mir]))
    return np.mean(agelst)

global ageham2fam

def corrcalc(sims):
    thegood = []
    agespread = []
    hammag = []
    for key in sims:
        if len(sims[key]) > 4:
            thegood.append(key)

    for key in thegood:
        for sublst in sims[key]:
            hammag.append(mag(sublst))
            agespread.append(avgage(fam2kids[key]))



    return agespread,hammag

def mag(lst):
    total = 0
    for el in lst:
        total+= el**2
    return max(lst)

regavg, regham = corrcalc(similarities)

print scs.spearmanr(regavg,regham)

plt.close()

plt.scatter(regham,regavg)
plt.show()


fle = open("textfiles/hsa_95fams_dollo_age-time.protein_list","r")
text = fle.readlines()
mirna2numage = {}
numage2howmany = {}
for line in text:
    if line[0] == "#":
        continue
    else:
        temp = string.replace(line,"\n","")
        temp = temp.split("\t")
        mirna2numage[temp[0]] = temp[1]

        if float(temp[1]) in numage2howmany.keys():
            numage2howmany[float(temp[1])]+=1
        else:
            numage2howmany[float(temp[1])] = 1


gaprate = [0]*(len(numage2howmany.keys()) - 1)

for i in range(0,len(numage2howmany.keys())-1):
    keys = numage2howmany.keys()
    keys.sort()
    gaprate[i] = float(numage2howmany[keys[i]]) / float(keys[i+1] - keys[i])

theages = numage2howmany.keys()
theages.sort()

for key in theages:
    key, "with", numage2howmany[key]



sortedbynum = sorted(numage2howmany.iteritems(), key=operator.itemgetter(1))


fle.close()
fle = open("textfiles/hsa_95fams_dollo_age-label.protein_list")

text = fle.readlines()
mirna2taxa = {}
taxa2num = {}
for line in text:
    if line[0] == "#":
        continue
    else:
        temp = string.replace(line,"\n","")
        temp = temp.split("\t")
        mirna2taxa[temp[0]] = temp[1]

        if temp[1] in taxa2num.keys():
            taxa2num[temp[1]]+=1
        else:
            taxa2num[temp[1]] = 1

newmiperyear = []
for i in range(0,len(theages)-1):
    deltaa = float(theages[i+1]) - float(theages[i])
    newmiperyear.append(float(numage2howmany[theages[i]]) / deltaa )

indexes = []
for i in range(0,len(newmiperyear)):
    indexes.append(i)
indexes.reverse()


subwords = {}

for di in diseases:
    temp = di.split(" ")
    for el in temp:
        if el in subwords:
            subwords[el] = subwords[el]+1
        else:
            subwords[el] = 1
bigclass = []
for key in subwords:
    if subwords[key] > 5 and key not in ["Disease", "Diseases", "Infection", "Disorders", "Chronic","Cell"]:
        bigclass.append(key)


fam2genvec = {}
spedis = bigclass
for di in diseases:
    temp = di.split(" ")
    good = True
    for el in temp:
        if el in bigclass:
            good = False

    if good == True:
        spedis.append(di)


spedisease2index = {}
for i in range(0,len(spedis)):
    spedisease2index[spedis[i]] = i


fam2kidspevec = {}
for fam in fam2kids:
    davecs = []
    for el in fam2kids[fam]:
        if el in mirna2disease:
            davecs.append(spebinvect(mirna2disease[el]))
    fam2kidspevec[fam] = davecs


spesims = {}


for key in fam2kidspevec:
    spesims[key] = vectable(fam2kidspevec[key])


spestd,speham = corrcalc(spesims)

print scs.stats.spearmanr(spestd,speham)


plt.close()
plt.scatter(speham,spestd)
plt.show()
plt.close()


def boxplotorg(mirna1):
    return avgage(fam2kids[mirna1])





simkeys  = similarities.keys()

simkeys.sort(key=boxplotorg)
badapples = []
for key in simkeys:
    if len(similarities[key]) < 4:
        badapples.append(key)

for key in badapples:
    simkeys.pop(simkeys.index(key))


boxplotitems = []

for key in simkeys:
    minilst = []
    for lst in similarities[key]:
        minilst.append( max(lst))
    boxplotitems.append(minilst)

plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)

ax.boxplot(boxplotitems)

fle.close()
fle = open("textfiles/fam2kids.txt","w")

for key in fam2kids:
    stringer = key
    for kid in fam2kids[key]:
        stringer+="\t"+kid
    fle.write(stringer+"\n")
fle.close
fle = open("textfiles/mirdisas.txt","w")
for key in mirna2disease:
    diseasestuff = mirna2disease[key]
    for thing in diseasestuff:
        fle.write(key+"\t"+thing+"\n")

print len(diseases)
print len(mirna2disease.keys())




