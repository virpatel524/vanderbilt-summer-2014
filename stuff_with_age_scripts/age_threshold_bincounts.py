# develop bincounts for the ages developed by PH for each consevation threshold using pyplot

__author__ = 'virpatel'
import stat


import numpy as np
from itertools import groupby
import math
from collections import Counter
lst = [50,60,75,80,85,95,100]
newlst = []
binlst = []
import matplotlib.pyplot as plt
ranger = range(0,17)
plots = []

txtname = ["hsa/hsa_","fams_dollo_age-depth.protein_list"]

for j in lst:
    binlst = []
    newlst = []
    text = open(txtname[0] + str(j) + txtname[1],"r")
    text = text.readlines()


    for i in text:
        if i[0]!= "#":
            temp = i.split(" ")
            temp = i.split("\t")
            temp2 = temp[1].split("\n")[0]
            temp = temp[0]
            newlst.append(temp2)

    for alpha in newlst:
        binlst.append(int(alpha))






    binlst.sort()

    freq = np.bincount(np.array(binlst),weights=None,minlength=17)





    pos = np.arange(len(ranger))

    width = 1.0
    ax = plt.axes()

    ax.set_xticks(pos + (width / 2))
    ax.set_xticklabels(ranger)

    plt.bar(pos,freq,width)
    plt.xlim(pos.min(),pos.max()+width)

    plt.savefig("bincounts/"+str(j)+"_binplot.png")

    plt.close()